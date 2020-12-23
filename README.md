# WideField Calcium Processing

The goal of this class and complementary functions is to process image files from a Wide-Field experiment and generate the fluorescence activity coming from it. 

## Download and Install

Clone `Widefield-calcium-post-processing` into your local directory and add it to the path. Also, in order to read the `.tif` files faster than Matlab's default library, `TIFFStack` repository is used. You will need to clone this one as well. The repository and further information can be found at:

[TIFFStack by Dilan Muir](https://github.com/DylanMuir/TIFFStack) 

## Usage

```
WF_obj = WideFieldProcessr('filename.tif')
```

The principal input of `WideFieldProcessor` is a mult-tif file with the extension `.tif`. The class will be called specifying the filename, present in the path, as input. If no input is selected, it will prompt the user to select one.


## Graphical methods






