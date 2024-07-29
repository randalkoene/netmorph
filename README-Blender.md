# Blender output in headless Blender mode

To run blender in headless mode in the background from Netmorph,
install any missing libraries:

```
sudo apt install libxxf86vm-dev
sudo apt install libxfixes3
sudo apt install libxi6
sudo apt install libxrender-dev libgl-dev
sudo apt install libxkbcommon-x11-0
sudo apt install libsm6
```

To add automatically producing OBJ and Blend output to any of
the example in the examples directory, add a second include
statement to the Netmorph call:

```
"include=../Source/examples/nesvbp/obj-and-blend"
```

To override the Blender executable path in this script, simply
add a "blender-exec-path=..." statement to the Netmorph call.

---
Randal A. Koene, 20240718
