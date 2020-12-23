# WideField Calcium Processing

The goal of this class and complementary functions is to process image files from a Wide-Field experiment and generate the fluorescence activity coming from it. 

## Download and Install

Clone `Widefield-calcium-post-processing` into your local directory and add it to the path. Also, in order to read the `.tif` files faster than Matlab's default library, `TIFFStack` repository is used. You will need to clone this one as well. The repository and further information can be found at:

[TIFFStack by Dilan Muir](https://github.com/DylanMuir/TIFFStack) 

## Usage

```
WF_obj = WideFieldProcessor('filename.tif')
```

The principal input of `WideFieldProcessor` is a mult-tif file with the extension `.tif`. In the case it is multiple tif files, function `TifsToStack.m` can be called to convert a series of single tifs into a multi-tif stack. The class will be called specifying the filename, present in the path, as input. If no input is specified, it will just create the `WF_obj` object, so that only the desired methods are applied. 

Additionally, to manually run the class with the filename specified:

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

### Specific frame reader

If the **Frame** option is set, then the class can read only the defined set of frames. Useful if the computer's RAM is not very large or if only certain chunks of the multi-tif are necessary, so that the processing time can be reduced. If the whole image is to be read, then the option `'Frame'` has to be set at `'off'` (default). Otherwise, the option value has to be a vector of two-dimensions.

```
WF_obj = WideFieldProcessor('filename.tif','Frame', 'off')
```

```
WF_obj.setOption('Frame', 1:100)
```

### Registration with reference image

If activated, this option lets register the average projection of the current multi-tiff with a reference image. If it exists, it will automatically look for the reference image of the current session `'[mouse_session]_refimg.tif'` and compare it to the mouse reference image `'[mouse]_reference.tif'`. If it can't find by name the reference of the session, it will prompt the user to choose a file or use the average projection instead. If it can't the mouse reference it will prompt the user to choose one.

```
WF_obj = WideFieldProcessor('filename.tif','RegisterSession', 'on')
```

```
WF_obj.setOption('RegisterSession', 'on')
```

```
[DFF_registered, tform] = WF_obj.registerSession(session_reference_filename, mouse_reference_filename)
```

#### Usage

A GUI will appear showing the session reference and the mouse reference, respectively. The user will then be prompted to choose the same points in both images (it has to be more than two). The transformation will then be applied to the whole DFF and a figure comparing the average projection and the reference will be showed. 

### Masking with parcellation template

If this option is chosed, then method `maskDFF(obj, template_struct)` will mask the DFF and remove the values outside the mask, so that the final DFF is lighter and can be saved faster. The input argument `template_struct` can be a structure or a matrix. In any case, the fields of the structure or the third dimension of the matrix have to be logical matrices of the same size as the individual frames of the stack. In order to create the template that matches a reference or the same stack, refere to  



## Graphical methods






