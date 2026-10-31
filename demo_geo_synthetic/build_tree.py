import hashlib
import json
import math

import numpy as np
from ete3 import Tree
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform

HERE = __file__.rsplit("/", 1)[0]
NE_JSON = HERE + "/../python_scripts/data/ne_110m_countries.json"

# Regional clusters of 3-5 countries. Within a cluster, countries are close
# on the tree (short, noisy branch distances) -- this is what makes 3-5
# countries chain together into a single connected group of map lines,
# instead of each pair forming its own isolated 2-country link.
#
# "eurasia_quad" also carries the Spain/Russia/Kazakhstan/Uzbekistan trio
# that must be pinned very close together regardless of noise (see below).
CLUSTERS = {
    "eurasia_quad": ["Spain", "Russia", "Kazakhstan", "Uzbekistan"],
    "west_europe": ["France", "Germany", "Italy", "Netherlands", "Belgium"],
    "iberian_atlantic": ["Portugal", "United Kingdom", "Morocco"],
    "north_europe": ["Sweden", "Norway", "Finland"],
    "east_europe": ["Ukraine", "Poland", "Romania"],
    "middle_east_ca": ["Turkey", "Iran", "Saudi Arabia", "Pakistan"],
    "east_asia": ["China", "Mongolia", "Japan", "South Korea", "Vietnam"],
    "southeast_asia": ["India", "Thailand", "Indonesia", "Philippines"],
    "north_america": ["United States of America", "Canada", "Mexico"],
    "central_america": ["Cuba", "Guatemala", "Panama", "Honduras", "Jamaica"],
    "south_america_andes": ["Peru", "Chile", "Bolivia", "Ecuador", "Colombia"],
    "south_america_atlantic": ["Brazil", "Argentina", "Venezuela", "Paraguay", "Uruguay"],
    "africa_west": ["Nigeria", "Ghana", "Senegal"],
    "africa_east": ["Kenya", "Tanzania", "Ethiopia"],
    "africa_north": ["Egypt", "Algeria", "Tunisia"],
    "oceania": ["Australia", "New Zealand", "Papua New Guinea"],
}

# Countries with samples but no close relative anywhere in the tree --
# each attaches on a long isolated branch, so they get a map dot but never
# a connection line, regardless of --top-k.
BACKGROUND_COUNTRIES = [
    "South Africa", "Angola", "Czechia", "Austria", "Switzerland", "Nepal",
    "Sri Lanka", "Malaysia", "Fiji", "Cambodia", "Zambia",
]

ALL_COUNTRIES = [c for members in CLUSTERS.values() for c in members] + BACKGROUND_COUNTRIES
assert len(ALL_COUNTRIES) == len(set(ALL_COUNTRIES)) == 72

# 1-3 samples per country (average 2), deterministic per-country (not
# random), so dots vary in size the way real sample sets would while
# keeping the total near the 2x density target (72 countries x ~2 = ~144).
SAMPLES_PER_COUNTRY = {c: 1 + (int(hashlib.sha256(c.encode()).hexdigest()[:4], 16) % 3)
                       for c in ALL_COUNTRIES}

# Scales haversine km -> tree branch-length units for the *cluster-level*
# UPGMA pass (cluster centroid vs. cluster centroid).
KM_TO_BRANCHLEN = 1.0 / 40000.0
TERMINAL_BRANCH = 0.003  # per-sample tip branch, with small deterministic jitter

# Within-cluster country distances: small, noisy, and close enough together
# that 3-5 countries chain into one connected component of map lines.
WITHIN_CLUSTER_BASE = 0.02
WITHIN_CLUSTER_NOISE_SIGMA = 0.35
WITHIN_CLUSTER_SEED = 1337

# Cluster-level (UPGMA) flatten + noise -- same purpose as before: avoid a
# perfectly clean "tight within continent / huge between continent" step
# function, and let an occasional cross-cluster pair sneak in.
FLATTEN_POWER = 0.55
NOISE_SIGMA = 0.22
RNG_SEED = 20261102

# Background countries attach this far from the rest of the tree -- far
# enough that they never qualify for a --top-k connection line.
BACKGROUND_BRANCH = 3.0

# Hand-tuned exact distances for the Spain/Russia/Kazakhstan/Uzbekistan
# quartet, grafted in at the eurasia_quad cluster's position in the
# noisy/flattened tree. Spain is much closer to Russia than to Kazakhstan,
# Uzbekistan sits a bit further out than Kazakhstan, and all of this is
# dramatically closer than any other within- or cross-cluster pair.
TRIO_SR_INTERNAL = 0.006    # spain<->russia common-ancestor depth
TRIO_KAZ_INTERNAL = 0.03    # kazakhstan's branch off the spain+russia clade
TRIO_UZB_INTERNAL = 0.05    # uzbekistan's branch off the spain+russia+kaz clade


def haversine_km(lat1, lon1, lat2, lon2):
    r = 6371.0
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlambda = math.radians(lon2 - lon1)
    a = math.sin(dphi / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dlambda / 2) ** 2
    return 2 * r * math.asin(math.sqrt(a))


def jitter(name, scale=0.0012):
    h = hashlib.sha256(name.encode()).hexdigest()
    frac = int(h[:8], 16) / 0xFFFFFFFF
    return scale * frac


def make_leaf_group(country, dist):
    group = Tree()
    group.dist = dist
    for k in range(1, SAMPLES_PER_COUNTRY[country] + 1):
        leaf_name = f"{country.replace(' ', '_')}_{k}"
        group.add_child(name=leaf_name, dist=TERMINAL_BRANCH + jitter(leaf_name))
    return group


def build_cluster_subtree(cluster_name, countries, rng):
    """UPGMA-free, deliberately shallow chain: sort countries by a
    deterministic-but-noisy key and nest them so adjacent countries in the
    chain are closest, giving a ladder where every adjacent pair -- and
    often several in a row -- lands inside any reasonable --top-k cutoff.
    This is what produces 3-5 country connected components instead of
    isolated pairs."""
    order = list(countries)
    # Deterministic shuffle per cluster so chain order isn't just input order.
    order.sort(key=lambda c: hashlib.sha256((cluster_name + c).encode()).hexdigest())

    node = make_leaf_group(order[0], 0.0)
    cumulative = 0.0
    for country in order[1:]:
        factor = math.exp(WITHIN_CLUSTER_NOISE_SIGMA * rng.standard_normal())
        step = WITHIN_CLUSTER_BASE * factor
        cumulative += step
        parent = Tree()
        parent.add_child(node, dist=step)
        parent.add_child(make_leaf_group(country, cumulative))
        node = parent
    return node, cumulative


def main():
    with open(NE_JSON) as fh:
        ne = json.load(fh)
    centroid = {c["name"]: (c["label_lat"], c["label_lon"]) for c in ne["countries"]}
    missing = [c for c in ALL_COUNTRIES if c not in centroid]
    if missing:
        raise SystemExit(f"missing countries in Natural Earth data: {missing}")

    cluster_names = list(CLUSTERS.keys())
    n = len(cluster_names)

    def cluster_centroid(members):
        lats = [centroid[c][0] for c in members]
        lons = [centroid[c][1] for c in members]
        return sum(lats) / len(lats), sum(lons) / len(lons)

    cluster_centroids = {name: cluster_centroid(members) for name, members in CLUSTERS.items()}

    km = np.zeros((n, n))
    for i in range(n):
        lat1, lon1 = cluster_centroids[cluster_names[i]]
        for j in range(i + 1, n):
            lat2, lon2 = cluster_centroids[cluster_names[j]]
            km[i, j] = km[j, i] = haversine_km(lat1, lon1, lat2, lon2)

    km_max = km.max()
    flat = km_max * np.power(np.clip(km / km_max, 0, 1), FLATTEN_POWER)

    rng = np.random.default_rng(RNG_SEED)
    noisy = flat.copy()
    for i in range(n):
        for j in range(i + 1, n):
            factor = math.exp(NOISE_SIGMA * rng.standard_normal())
            noisy[i, j] = noisy[j, i] = flat[i, j] * factor

    dist = noisy * KM_TO_BRANCHLEN

    condensed = squareform(dist, checks=False)
    Z = linkage(condensed, method="average")

    within_rng = np.random.default_rng(WITHIN_CLUSTER_SEED)
    cluster_subtree = {}
    cluster_depth = {}
    for name in cluster_names:
        subtree, depth = build_cluster_subtree(name, CLUSTERS[name], within_rng)
        cluster_subtree[name] = subtree
        cluster_depth[name] = depth

    cluster_height = dict(enumerate([0.0] * n))
    cluster_node = {}
    for row_i, (a, b, height, _count) in enumerate(Z):
        a, b = int(a), int(b)
        new_id = n + row_i
        internal = Tree()
        for child_id in (a, b):
            child_height = cluster_height.get(child_id, 0.0)
            if child_id < n:
                name = cluster_names[child_id]
                branch_len = max(height - cluster_depth[name], 0.0)
                internal.add_child(cluster_subtree[name], dist=branch_len)
            else:
                branch_len = max(height - child_height, 0.0)
                internal.add_child(cluster_node[child_id], dist=branch_len)
        cluster_node[new_id] = internal
        cluster_height[new_id] = height

    final_id = n + len(Z) - 1
    root = cluster_node[final_id]

    # --- graft Spain/Russia/Kazakhstan/Uzbekistan onto the eurasia_quad
    # cluster's natural position, overriding the generic within-cluster
    # chain built above for just these four countries. ---------------------
    eurasia_subtree = cluster_subtree["eurasia_quad"]
    eurasia_parent = eurasia_subtree.up
    eurasia_attach_dist = eurasia_subtree.dist
    eurasia_subtree.detach()

    spain_group = make_leaf_group("Spain", TRIO_SR_INTERNAL)
    russia_group = make_leaf_group("Russia", TRIO_SR_INTERNAL)
    sr = Tree()
    sr.add_child(spain_group)
    sr.add_child(russia_group)

    kaz_group = make_leaf_group("Kazakhstan", TRIO_KAZ_INTERNAL)
    srk = Tree()
    srk.add_child(sr, dist=TRIO_KAZ_INTERNAL - TRIO_SR_INTERNAL)
    srk.add_child(kaz_group)

    uzb_group = make_leaf_group("Uzbekistan", TRIO_UZB_INTERNAL)
    quad = Tree()
    quad.dist = eurasia_attach_dist
    quad.add_child(srk, dist=TRIO_UZB_INTERNAL - TRIO_KAZ_INTERNAL)
    quad.add_child(uzb_group)

    eurasia_parent.add_child(quad)

    # --- attach background countries on long, isolated branches ----------
    for country in BACKGROUND_COUNTRIES:
        root.add_child(make_leaf_group(country, BACKGROUND_BRANCH))

    out_nwk = HERE + "/synthetic_countries.nwk"
    root.write(outfile=out_nwk, format=1)

    # --- write matching sample_locations.tsv ------------------------------
    tsv_path = HERE + "/sample_locations.tsv"
    lines = ["genome_id\tlocation"]
    for country in ALL_COUNTRIES:
        stem = country.replace(" ", "_")
        for k in range(1, SAMPLES_PER_COUNTRY[country] + 1):
            lines.append(f"{stem}_{k}\t{country}")
    with open(tsv_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")

    total_samples = sum(SAMPLES_PER_COUNTRY.values())
    print(f"{len(ALL_COUNTRIES)} countries, {total_samples} samples")

    checks = [
        ("Spain_1", "Russia_1", "Spain-Russia (grafted, tightest)"),
        ("Spain_1", "Kazakhstan_1", "Spain-Kazakhstan (grafted, close)"),
        ("Spain_1", "Uzbekistan_1", "Spain-Uzbekistan (grafted, close-ish)"),
        ("Russia_1", "Kazakhstan_1", "Russia-Kazakhstan (grafted)"),
        ("France_1", "Germany_1", "west_europe intra-cluster"),
        ("France_1", "Italy_1", "west_europe intra-cluster"),
        ("Brazil_1", "Argentina_1", "south_america_atlantic intra-cluster"),
        ("Cuba_1", "Panama_1", "central_america intra-cluster"),
        ("Nigeria_1", "Ghana_1", "africa_west intra-cluster"),
        ("South_Africa_1", "Zambia_1", "background-background (should be far)"),
        ("Spain_1", "China_1", "cross-cluster (far)"),
        ("Nigeria_1", "Brazil_1", "cross-cluster (far)"),
    ]
    print()
    for a, b, label in checks:
        try:
            print(f"{label:45s} {a:30s} {b:30s} {root.get_distance(a, b):.4f}")
        except Exception as e:
            print(label, a, b, "ERROR", e)


if __name__ == "__main__":
    main()
