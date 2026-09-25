from __future__ import annotations

import os
from pathlib import Path
import shutil
import subprocess

from Cython.Build import cythonize
from setuptools import Extension, find_packages, setup
from setuptools.command.build_ext import build_ext
from setuptools.command.build_py import build_py


root = Path(__file__).resolve().parent
version = os.environ.get("POLYPAVES_VERSION", "0.3.0")


class StrippedBuildExt(build_ext):
    @staticmethod
    def _portable_command(command):
        cleaned = []
        skip_next = False
        for argument in command:
            if skip_next:
                skip_next = False
                continue
            if argument == "-B":
                skip_next = True
                continue
            if argument.startswith(("-L/home/", "-Wl,-rpath,/home/", "-Wl,-rpath-link,/home/")):
                continue
            cleaned.append(argument)
        return cleaned

    def run(self):
        super().run()
        strip = shutil.which("strip")
        if strip:
            for extension in self.extensions:
                subprocess.run(
                    [strip, "--strip-unneeded", self.get_ext_fullpath(extension.name)],
                    check=True,
                )

    def build_extensions(self):
        self.compiler.compiler_so = self._portable_command(self.compiler.compiler_so)
        self.compiler.linker_so = self._portable_command(self.compiler.linker_so)
        super().build_extensions()


# Every module of the package except the facade is compiled; none ships as source.
COMPILED = sorted(path.stem for path in (root / "polypaves").glob("*.py") if path.stem != "__init__")


class PublicBuildPy(build_py):
    def find_package_modules(self, package, package_dir):
        modules = super().find_package_modules(package, package_dir)
        return [module for module in modules if module[1] not in COMPILED]


extensions = cythonize(
    [
        Extension(
            f"polypaves.{name}",
            [f"polypaves/{name}.py"],
            extra_compile_args=[
                "-O3",
                "-fvisibility=hidden",
                f"-ffile-prefix-map={root}=.",
                f"-fdebug-prefix-map={root}=.",
            ],
        )
        for name in COMPILED
    ],
    compiler_directives={
        "binding": True,
        "embedsignature": False,
        "emit_code_comments": False,
        "language_level": 3,
        "linetrace": False,
        "profile": False,
    },
)

setup(
    name="polypaves",
    version=version,
    description="PAVES: Polymer Automation for Virtual Evaluation and Simulation",
    python_requires=">=3.10",
    packages=find_packages(),
    package_data={"polypaves": ["_native*.so"]},
    include_package_data=False,
    ext_modules=extensions,
    cmdclass={"build_ext": StrippedBuildExt, "build_py": PublicBuildPy},
    zip_safe=False,
    # `polypaves-build INPUT.paves [key=value ...]`: the file-driven builder, from the compiled cli module.
    entry_points={"console_scripts": ["polypaves-build = polypaves.cli:main"]},
)
