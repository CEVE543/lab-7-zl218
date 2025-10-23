#!/usr/bin/env julia
"""
Download CMIP6 data for Boston from Google Cloud Storage and save as NetCDF4
"""

using Pkg
Pkg.activate(@__DIR__)
using YAXArrays, Zarr, DimensionalData, NetCDF

# Boston coordinates
const BOSTON_LAT = 42.3631
const BOSTON_LON = 360 - 71.0064  # Convert to 0-360

# CMIP6 URLs - using GFDL-ESM4 which has both historical and future runs
const HIST_URL = "gs://cmip6/CMIP6/CMIP/NOAA-GFDL/GFDL-ESM4/historical/r1i1p1f1/3hr/tas/gr1/v20190726"
const SSP370_URL = "gs://cmip6/CMIP6/ScenarioMIP/NOAA-GFDL/GFDL-ESM4/ssp370/r1i1p1f1/3hr/tas/gr1/v20180701"

"""Download data from GCS and save to NetCDF4"""
function download_and_save(url::String, output_file::String, lat::Float64, lon::Float64)
    # Open remote store and get temperature data
    ds = open_dataset(zopen(url, consolidated=true))
    tas = ds["tas"]

    # Select nearest grid cell (this still requires downloading spatial chunks)
    println("      Selecting nearest grid cell...")
    boston_data = tas[lon=At(lon, atol=2.0), lat=At(lat, atol=2.0)]

    # Save directly to NetCDF using YAXArrays
    # The download happens during the save operation
    println("      Downloading and saving to NetCDF...")
    savecube(boston_data, output_file, driver=:netcdf)

    return filesize(output_file) / 1024^2  # Return size in MB
end

function main()
    println("=" ^ 80)
    println("Downloading CMIP6 data for Boston Logan Airport")
    println("Location: $(BOSTON_LAT)°N, $(BOSTON_LON)°E")
    println("Model: NOAA-GFDL GFDL-ESM4")
    println("Variable: Near-surface air temperature (tas)")
    println("Temporal resolution: 3-hourly")
    println("=" ^ 80)

    # Download historical
    hist_file = joinpath(@__DIR__, "boston_historical.nc")
    println("\n1. Downloading historical scenario (1850-2014)...")
    println("   Ensemble member: r1i1p1f1")
    println("   Progress: ")
    flush(stdout)
    hist_size = download_and_save(HIST_URL, hist_file, BOSTON_LAT, BOSTON_LON)
    println("✓ $(round(hist_size, digits=1)) MB")

    # Download SSP3-7.0
    ssp370_file = joinpath(@__DIR__, "boston_ssp370.nc")
    println("\n2. Downloading SSP3-7.0 scenario (2015-2100)...")
    println("   Ensemble member: r1i1p1f1")
    println("   Progress: ")
    flush(stdout)
    ssp370_size = download_and_save(SSP370_URL, ssp370_file, BOSTON_LAT, BOSTON_LON)
    println("✓ $(round(ssp370_size, digits=1)) MB")

end

main()
