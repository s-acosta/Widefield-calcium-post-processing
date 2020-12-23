# WideField Calcium Processing

The goal of this class and complementary functions is to process image files from a Wide-Field experiment and generate the fluorescence activity coming from it. 

## Download and Install

Clone `Widefield-calcium-post-processing` into your local directory and add it to the path. Also, in order to read the `.tif` files faster than Matlab's default library, `TIFFStack` repository is used. You will need to clone this one as well. The repository and further information can be found at:

[TIFFStack by Dilan Muir](https://github.com/DylanMuir/TIFFStack) 

## Usage

```
WF_obj = WideFieldProcessor('filename.tif')
```

The principal input of `WideFieldProcessor` is a mult-tif file with the extension `.tif`. The class will be called specifying the filename, present in the path, as input. If there is no input argument, it will prompt the user to select one. 

To run the class manually:

```
WF_obj = WideFieldProcessor('filename.tif','Mode','manual')
```


## Parameters

The class accepts two ways of accepting change of options. If it is run in auto mode (default mode):

```
WF_obj = WideFieldProcessor('filename.tif','OptionName', optionValue)
```

The defaults names and values of the parameters can be found below. In the manual mode (and also in the auto mode if we want to change anything after initializing), parameters can be changed using:

```
WF_obj.setOption('OptionName', optionValue)
```

## Features




## Graphical methods






