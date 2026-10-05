# Split Singularity build into base + app images (deferred)

## Problem

`singularity build` has no layer cache for `%post`/`%files`/`%setup` — every build
re-executes the ENTIRE `%post` script from scratch against a fresh copy of the base
rootfs (`Bootstrap`/`From:`), no matter how small the change. Only the base image
pull itself is cached (visible via `singularity cache list` as `library`/`blob`
entries for `mambaorg/micromamba`).

For `Singularity.OP.v3.1.0.sing_test.recipe`, this means every rebuild re-runs, in
full, every time, even for a one-line change (e.g. the `python=3.12` pin added
2026-09-28 for the ete3/`cgi`-removal fix):
- the `micromamba install`/`micromamba create` dependency solves (slow)
- the CheckM2 DIAMOND DB download (`checkm2 database --download`, multi-GB)
- the ASTRAL / catfasta2phyml / Alignment_Assessment git clones

This is a real architectural limitation of the native SIF build model (Sylabs
prioritizes reproducible flat images over incremental build speed), not a config
bug — checked `/etc/singularity/singularity.conf`, no caching knob changes it.

## Fix plan: base image + app image

Split the recipe into two definition files:

1. **`Singularity.OP.base.recipe`** — everything expensive and rarely-changing:
   apt packages, micromamba installs (both `base` and `gather_genomes` envs),
   CheckM2 DB download, ASTRAL/catfasta2phyml/Alignment_Assessment clones.
   Only rebuild this when a dependency/version actually changes.

2. **`Singularity.OP.app.recipe`** — bootstraps FROM the base image instead of
   `docker://mambaorg/micromamba`:
   ```
   Bootstrap: localimage
   From: OrthoPhyl_base.sif
   ```
   Only does the fast part: clone the OrthoPhyl repo/branch, `chmod`,
   `%environment`/`%runscript`/`%help`. Rebuilding this after a code/branch change
   is seconds, not the full multi-minute dependency-solve + DB-download cycle.

## Caveat specific to this host

Local `--fakeroot` builds are blocked here — `/proc/sys/user/max_user_namespaces`
is `0` (user namespaces disabled kernel-wide; a sysadmin would need to
`sysctl -w user.max_user_namespaces=15000` and persist it under `/etc/sysctl.d/` to
unblock this). So builds must go through `--remote` (Sylabs Cloud remote builder,
already configured — `singularity remote list` shows `SylabsCloud`).

`Bootstrap: localimage` won't see a file that only exists on your machine when
building remotely. Workaround:
1. Build the base image once via `--remote`.
2. `singularity push` it to your personal Sylabs Cloud library
   (e.g. `library://<you>/default/orthophyl-base:latest`).
3. Point the app recipe's `Bootstrap`/`From:` at that library reference instead of
   `localimage`.

Then only the app-stage recipe needs a remote-build round-trip for routine
branch/permission changes; the expensive base stays untouched until its own
recipe changes.

## Status

Not started. Revisit after the current CheckM2-channel-solve build error
(`checkm2 =* * does not exist`) is resolved and the immediate `python=3.12`-pinned
build succeeds.
