# 2022_Cimini_NatureProtocols

## Download cell images

Cell images are available on an S3 bucket.
The images can be downloaded using the command

```bash
plate = BR00117011__2020-11-08T19_57_47-Measurement1
aws s3 cp \
  --no-sign-request \
  --recursive \
  s3://cellpainting-gallery/jump-pilot/source_4/images/2020_11_04_CPJUMP1/images/${plate} .
```

You can test out download for a single file using:

```
suffix=BR00117011__2020-11-08T19_57_47-Measurement1/Images/r01c01f01p01-ch1sk1fk1fl1.tiff

aws s3 cp \
  --no-sign-request \
  s3://cellpainting-gallery/jump-pilot/source_4/images/2020_11_04_CPJUMP1/images/${suffix} .
```

## Download illumination correction files

Illumination correction files, produced by running illum.cppipe on the cell images, are available on an S3 bucket. 
The images can be downloaded using the command

```bash
plate = BR00117011
aws s3 cp \
  --no-sign-request \
  --recursive \
  s3://cellpainting-gallery/jump-pilot/source_4/images/2020_11_04_CPJUMP1/illum/${plate} .
```

## Image metadata

The image metadata will be automatically extracted from the images using the pipelines provided.
An orientation to the metadata is as follows:

The folder for the 384-well plate contains images from nine sites for each well.
There are eight images per site (five from the fluorescent channels and three brightfield images).
The names of the image files follow the naming convention - `rXXcXXfXXp01-chXXsk1fk1fl1.tiff` where
- `rXX` is the row number of the well that was imaged. `rXX` ranges from `r01` to `r16`.
- `cXX` is the column number of the well that was imaged. `cXX` ranges from `c01` to `c24`.
- `fXX` corresponds to the site that was imaged. `fXX` ranges from `f01` to `f09`.
- `chXX` corresponds to the fluorescent channels imaged. `chXX` ranges from `ch01` to `ch08`.
    - `ch01` - Alexa 647
    - `ch02` - Alexa 568
    - `ch03` - Alexa 488 long
    - `ch04` - Alexa 488
    - `ch05` - Hoechst 33342
    - `ch06-8` - three brightfield z planes.
