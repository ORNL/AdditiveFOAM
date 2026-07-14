## Nix

First install the [Nix package manager][NIX] and then enable
[Flakes][Flakes].  Alternatively, check out the [Determinate Systems
Installer][Determinate] for an out of the box experience. See
[nix.dev][nix.dev] for more help with Nix.

To get an environment with both AdditiveFOAM and ExaCA use

    $ nix develop github:ORNL/AdditiveFOAM?dir=envs/nix --accept-flake-config
    # AdditiveFOAM and ExaCA now available
    $ additiveFoam --help
	$ ExaCA --help

To get a specific release, use:

    $ nix develop github:ORNL/AdditiveFOAM/<tag>?dir=envs/nix --accept-flake-config

In principle the build should be fast (on x86-64-linux) as the
binaries are downloaded from Cachix. You will be prompted to `allow
configuration setting 'extra-substituters'`.

See the main [README](../../README.md) for more details on building
AdditiveFOAM.

[NIX]: https://nixos.org/download.html
[Flakes]: https://nixos.wiki/wiki/Flakes
[nix.dev]: https://nix.dev
[Determinate]: https://github.com/DeterminateSystems/nix-installer

