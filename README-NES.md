# Using Netmorph as an embedded component of NES

As of 2024-07-06, Netmorph is embedded for use within the Neural Emulation System (NES)
and Virtual Brain Platform (VBP).

This depends on the following:

1. An alternative entry point through `Include/Netmorph.cpp:main2()` (instead of `nibr.cc:main()`).
2. Redirection of the logging methods `error()`, `warning()`, `report()`, `progress()` in `global.hh`
   and a recommendation to avoid any direct use of `std::cout` or `std::cerr`.
3. Netmorph closing post-op that does not exit all programs.
4. Netmorph starting on a thread of an NES API handler, using the `std::string` method of supplying
   model commands to Netmorph through the `Command_Line_Parameters` object.

rkoene@carboncopies.org
