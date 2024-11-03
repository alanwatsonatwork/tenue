import tenue.instrument

tenue.instrument.setvalues(
    overscanyslice=slice(1, 10),
    overscanxslice=slice(18, 2065),
    trimyslice=slice(11, 2058),
    trimxslice=slice(18, 2065),
    filterkeyword="FILTER",
    alphakeyword="SMTMRA",
    deltakeyword="SMTMDE",
    pixelscale=0.40 / 3600,
)
