## Nix

First install the [Nix package manager][NIX] and then enable
[Flakes][Flakes].  Alternatively, check out the [Determinate Systems
Installer][Determinate] for an out of the box experience. See
[nix.dev][nix.dev] for more help with Nix.

To get an environment with both AdditiveFOAM and ExaCA use

    $ nix develop github:ORNL/AdditiveFOAM
    # AdditiveFOAM and ExaCA now available
    $ additiveFoam --help
	$ ExaCA --help

[NIX]: https://nixos.org/download.html
[Flakes]: https://nixos.wiki/wiki/Flakes
[nix.dev]: https://nix.dev
[Determinate]: https://github.com/DeterminateSystems/nix-installer
