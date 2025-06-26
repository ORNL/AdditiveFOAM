# Taken from https://github.com/NixOS/nixpkgs/blob/5ee46705af97c806817aa41dbb2914a4bf431ee3/pkgs/by-name/op/openfoam-com/generic.nix

{
  flex,
  bison,
  scotch,
  zlib,
  fetchFromGitLab,
  mpi
}:
scotch.overrideAttrs (prev: rec {
	version = "6.0.9";

  src = fetchFromGitLab {
    domain = "gitlab.inria.fr";
    owner = "scotch";
    repo = "scotch";
    rev = "v${version}";
    hash = "sha256-YR2Evj7VnDNGQVNlHovj6CfvPJgAoHeTKu4EW5cXtdg=";
  };

  nativeBuildInputs = [ ];

  buildInputs = [
    bison
    mpi
    flex
    zlib
  ];

  preConfigure = ''
    cd src
    ln -s Make.inc/Makefile.inc.x86-64_pc_linux2 Makefile.inc
  '';

  buildFlags = [ "scotch ptscotch" ];

  installFlags = [ "prefix=\${out}" ];
})
