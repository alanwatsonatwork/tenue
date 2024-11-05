import tenue.instrument

tenue.instrument.setvalues(
    datamax=65535,
    overscanyslice=slice(1, 10),
    overscanxslice=slice(18, 2065),
    trimyslice=slice(11, 2058),
    trimxslice=slice(18, 2065),
    filterkeyword="FILTER",
    exposuretimekeyword="EXPTIME",
    alphakeyword="SMTMRA",
    deltakeyword="SMTMDE",
    rotationkeyword="SMTMRO",
    pixelscale=0.40 / 3600,
    rotation=63.703,
    flatmax=32000,
)
