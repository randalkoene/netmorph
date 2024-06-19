# NESVBP

This example builds a model with 7 neurons in 3 layers, 2 of which are
interneurons. This is meant to test building a model architecture
similar to the XOR example in VBP.

To run the simulation, execut the `run.sh` script.

To include the resulting model in a VBP/NES model, copy the resulting
output files to a folder accessible by VBP (if necessary). Then,
in VBP, run the ground-truth model script in the xor_scnm model
directory of VBP and use the `-NmSource` argument to point it to
the common trunk of the output files, for example:

```
./xor_scnm_groundtruth.py -NmSource ./netmorph_output_data/nesvbp/nesvbp_20240618
```

---
Randal A. Koene, 20240618
