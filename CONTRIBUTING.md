# Contributing

Contributing to AdditiveFOAM is easy: just open a
[pull request](https://help.github.com/articles/using-pull-requests/).
Make `main` the destination branch on the [AdditiveFOAM
repository](https://github.com/ORNL/AdditiveFOAM) and allow edits from
maintainers.

Your pull request must work with all current AdditiveFOAM tutorial examples
and be reviewed by at least one AdditiveFOAM developer.

For local verification, build the code with `./Allwmake`, build the native
test harness with `./tests/Allwmake`, and run it with `./tests/run`. The test
workflow and instructions for adding new native tests are documented in
[TESTING.md](TESTING.md).
