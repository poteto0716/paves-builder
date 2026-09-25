# Third-party materials

polypaves binaries, documentation, examples and runtime scripts are distributed under
[LICENSE](LICENSE). This does not grant rights to third-party materials.

## Compiler and system runtime

The executable was compiled with the standard Ubuntu GCC 13.3 toolchain. It uses
libstdc++, libgcc, glibc/libm and the ELF loader dynamically. No third-party shared
library or static library is shipped separately. CMake links only the project's
own core library into the executable; no EMC executable or recovered library is
built or linked. Internal probe/test applications are not distributed.

GCC runtime/header portions are subject to GPLv3 with the GCC Runtime Library
Exception 3.1. The exception covers eligible compilation of independent modules,
including proprietary programs. See the
[official libstdc++ license documentation](https://gcc.gnu.org/onlinedocs/libstdc++/manual/license.html),
[GPLv3](LICENSES/GPL-3.0.txt) and
[exception](LICENSES/GCC-Runtime-Library-Exception-3.1.txt).
System glibc is installed by Ubuntu, not redistributed here; its notices and
terms are available in the system package documentation. Python, Bash and LAMMPS
are also user-installed external dependencies, not bundled software.

## PCFF / IFF / silica

These databases and structure files are included under `external/` so the examples
are directly reproducible. The inspected EMC installation contains GPLv3 text;
the official project's listing says MIT.
The historical PCFF database has its own scientific provenance, and its scope of
redistribution permission was not established by either indication.

The IFF 1.5 archive and the authors' public GitHub repository contain parameters
and surface models, but the inspected material does not state an explicit
redistribution license. Citation requests are recorded; they are not treated as
a change to the files' terms. The bundled files remain third-party materials and
are not licensed under the polypaves Binary Software License.

[External data](docs/external_data.md) gives official sources, exact filenames,
installation locations, hashes and citations. The installer can restore the
pinned copies from their original distributions. Generated LAMMPS files remain
excluded from version control.

## Scope of this audit

Build configuration, dependencies, the included runtime helpers and data provenance
were inspected. No license is inferred from a file merely being present in the
development installation. No source-available or open-source license is granted
to the polypaves core. The copyright holder supplied the license text in this release.
