## See NIX.md for help getting started with Nix

{
  description = "An open-source CFD code for additive manufacturing built on OpenFOAM.";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixos-25.05";
    exaca.url   = "github:LLNL/ExaCA?dir=envs/nix";
    parts.url = "github:hercules-ci/flake-parts";
  };

  outputs = inputs @ { self, parts, ... }: (
    parts.lib.mkFlake { inherit inputs; } {
      systems = [
        "x86_64-linux"
      ];

      perSystem = { pkgs, inputs', ... }: {
        packages = rec {
          default = additivefoam;
            
          openfoam = pkgs.callPackage ./nix/openfoam.nix { scotch = scotch; };
          scotch = pkgs.callPackage ./nix/scotch.nix { };
          exaca = inputs'.exaca.packages.default;

          additivefoam = pkgs.callPackage ./nix/additivefoam.nix {
            src = self;
            version = "master";
            exaca = exaca;
            openfoam = openfoam;
          };
        };

        imports = [
          ./nix/dev.nix
        ];
      };
    }
  );
}
