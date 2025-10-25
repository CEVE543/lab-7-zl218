// Some definitions presupposed by pandoc's typst output.
#let blockquote(body) = [
  #set text( size: 0.92em )
  #block(inset: (left: 1.5em, top: 0.2em, bottom: 0.2em))[#body]
]

#let horizontalrule = line(start: (25%,0%), end: (75%,0%))

#let endnote(num, contents) = [
  #stack(dir: ltr, spacing: 3pt, super[#num], contents)
]

#show terms: it => {
  it.children
    .map(child => [
      #strong[#child.term]
      #block(inset: (left: 1.5em, top: -0.4em))[#child.description]
      ])
    .join()
}

// Some quarto-specific definitions.

#show raw.where(block: true): set block(
    fill: luma(230),
    width: 100%,
    inset: 8pt,
    radius: 2pt
  )

#let block_with_new_content(old_block, new_content) = {
  let d = (:)
  let fields = old_block.fields()
  fields.remove("body")
  if fields.at("below", default: none) != none {
    // TODO: this is a hack because below is a "synthesized element"
    // according to the experts in the typst discord...
    fields.below = fields.below.abs
  }
  return block.with(..fields)(new_content)
}

#let empty(v) = {
  if type(v) == str {
    // two dollar signs here because we're technically inside
    // a Pandoc template :grimace:
    v.matches(regex("^\\s*$")).at(0, default: none) != none
  } else if type(v) == content {
    if v.at("text", default: none) != none {
      return empty(v.text)
    }
    for child in v.at("children", default: ()) {
      if not empty(child) {
        return false
      }
    }
    return true
  }

}

// Subfloats
// This is a technique that we adapted from https://github.com/tingerrr/subpar/
#let quartosubfloatcounter = counter("quartosubfloatcounter")

#let quarto_super(
  kind: str,
  caption: none,
  label: none,
  supplement: str,
  position: none,
  subrefnumbering: "1a",
  subcapnumbering: "(a)",
  body,
) = {
  context {
    let figcounter = counter(figure.where(kind: kind))
    let n-super = figcounter.get().first() + 1
    set figure.caption(position: position)
    [#figure(
      kind: kind,
      supplement: supplement,
      caption: caption,
      {
        show figure.where(kind: kind): set figure(numbering: _ => numbering(subrefnumbering, n-super, quartosubfloatcounter.get().first() + 1))
        show figure.where(kind: kind): set figure.caption(position: position)

        show figure: it => {
          let num = numbering(subcapnumbering, n-super, quartosubfloatcounter.get().first() + 1)
          show figure.caption: it => {
            num.slice(2) // I don't understand why the numbering contains output that it really shouldn't, but this fixes it shrug?
            [ ]
            it.body
          }

          quartosubfloatcounter.step()
          it
          counter(figure.where(kind: it.kind)).update(n => n - 1)
        }

        quartosubfloatcounter.update(0)
        body
      }
    )#label]
  }
}

// callout rendering
// this is a figure show rule because callouts are crossreferenceable
#show figure: it => {
  if type(it.kind) != str {
    return it
  }
  let kind_match = it.kind.matches(regex("^quarto-callout-(.*)")).at(0, default: none)
  if kind_match == none {
    return it
  }
  let kind = kind_match.captures.at(0, default: "other")
  kind = upper(kind.first()) + kind.slice(1)
  // now we pull apart the callout and reassemble it with the crossref name and counter

  // when we cleanup pandoc's emitted code to avoid spaces this will have to change
  let old_callout = it.body.children.at(1).body.children.at(1)
  let old_title_block = old_callout.body.children.at(0)
  let old_title = old_title_block.body.body.children.at(2)

  // TODO use custom separator if available
  let new_title = if empty(old_title) {
    [#kind #it.counter.display()]
  } else {
    [#kind #it.counter.display(): #old_title]
  }

  let new_title_block = block_with_new_content(
    old_title_block, 
    block_with_new_content(
      old_title_block.body, 
      old_title_block.body.body.children.at(0) +
      old_title_block.body.body.children.at(1) +
      new_title))

  block_with_new_content(old_callout,
    block(below: 0pt, new_title_block) +
    old_callout.body.children.at(1))
}

// 2023-10-09: #fa-icon("fa-info") is not working, so we'll eval "#fa-info()" instead
#let callout(body: [], title: "Callout", background_color: rgb("#dddddd"), icon: none, icon_color: black, body_background_color: white) = {
  block(
    breakable: false, 
    fill: background_color, 
    stroke: (paint: icon_color, thickness: 0.5pt, cap: "round"), 
    width: 100%, 
    radius: 2pt,
    block(
      inset: 1pt,
      width: 100%, 
      below: 0pt, 
      block(
        fill: background_color, 
        width: 100%, 
        inset: 8pt)[#text(icon_color, weight: 900)[#icon] #title]) +
      if(body != []){
        block(
          inset: 1pt, 
          width: 100%, 
          block(fill: body_background_color, width: 100%, inset: 8pt, body))
      }
    )
}



#let article(
  title: none,
  subtitle: none,
  authors: none,
  date: none,
  abstract: none,
  abstract-title: none,
  cols: 1,
  lang: "en",
  region: "US",
  font: "libertinus serif",
  fontsize: 11pt,
  title-size: 1.5em,
  subtitle-size: 1.25em,
  heading-family: "libertinus serif",
  heading-weight: "bold",
  heading-style: "normal",
  heading-color: black,
  heading-line-height: 0.65em,
  sectionnumbering: none,
  toc: false,
  toc_title: none,
  toc_depth: none,
  toc_indent: 1.5em,
  doc,
) = {
  set par(justify: true)
  set text(lang: lang,
           region: region,
           font: font,
           size: fontsize)
  set heading(numbering: sectionnumbering)
  if title != none {
    align(center)[#block(inset: 2em)[
      #set par(leading: heading-line-height)
      #if (heading-family != none or heading-weight != "bold" or heading-style != "normal"
           or heading-color != black) {
        set text(font: heading-family, weight: heading-weight, style: heading-style, fill: heading-color)
        text(size: title-size)[#title]
        if subtitle != none {
          parbreak()
          text(size: subtitle-size)[#subtitle]
        }
      } else {
        text(weight: "bold", size: title-size)[#title]
        if subtitle != none {
          parbreak()
          text(weight: "bold", size: subtitle-size)[#subtitle]
        }
      }
    ]]
  }

  if authors != none {
    let count = authors.len()
    let ncols = calc.min(count, 3)
    grid(
      columns: (1fr,) * ncols,
      row-gutter: 1.5em,
      ..authors.map(author =>
          align(center)[
            #author.name \
            #author.affiliation \
            #author.email
          ]
      )
    )
  }

  if date != none {
    align(center)[#block(inset: 1em)[
      #date
    ]]
  }

  if abstract != none {
    block(inset: 2em)[
    #text(weight: "semibold")[#abstract-title] #h(1em) #abstract
    ]
  }

  if toc {
    let title = if toc_title == none {
      auto
    } else {
      toc_title
    }
    block(above: 0em, below: 2em)[
    #outline(
      title: toc_title,
      depth: toc_depth,
      indent: toc_indent
    );
    ]
  }

  if cols == 1 {
    doc
  } else {
    columns(cols, doc)
  }
}

#set table(
  inset: 6pt,
  stroke: none
)
#import "@preview/fontawesome:0.5.0": *

#set page(
  paper: "us-letter",
  margin: (x: 1in,y: 1in,),
  numbering: "1",
)

#show: doc => article(
  title: [CEVE 543 Fall 2025 Lab 7: Bias Correction Implementation],
  subtitle: [Delta method and quantile mapping for temperature bias correction],
  authors: (
    ( name: [Zijie (Ian) Liang],
      affiliation: [],
      email: [] ),
    ),
  date: [2025-10-24],
  fontsize: 11pt,
  sectionnumbering: "1.1.a",
  toc_title: [Table of contents],
  toc_depth: 3,
  cols: 1,
  doc,
)

= Background and Goals
<background-and-goals>
Climate models have systematic biases that directly affect impact assessments, arising from coarse spatial resolution, parameterization of sub-grid processes, representation of topography, and errors in simulated circulation patterns. This lab implements two widely-used bias correction approaches: the delta method and quantile-quantile (QQ) mapping. The delta method preserves the climate model's change signal while anchoring absolute values to observations, whereas QQ-mapping corrects the full distribution of values.

Both methods assume stationarity---that the statistical relationship between model and observations remains constant across climate states. This assumption may not hold under significant climate change. We'll explore the strengths and limitations of each method using temperature data for Boston, providing hands-on experience before PS2 Part 1.

= Study Location and Data
<study-location-and-data>
This lab uses temperature data for Boston Logan International Airport (Station ID: USW00014739, 42.3631$""^circle.stroked.tiny$N, 71.0064$""^circle.stroked.tiny$W). Observational data comes from GHCN-Daily (1936-2024) and is provided in the `USW00014739.csv` file in degrees Celsius. Climate model data comes from the GFDL-ESM4 model's 3-hourly near-surface air temperature (`tas`), pre-downloaded from Google Cloud Storage for both historical (1850-2014) and SSP3-7.0 (2015-2100) scenarios. Refer to #link(".\\labs/Lab-6/index.qmd")[Lab 6] for details on CMIP6 data structure.

== Data Processing Notes
<data-processing-notes>
The GHCN data provides daily average temperature in degrees Celsius. CMIP6 provides 3-hourly instantaneous temperature in Kelvin, which we'll convert to Celsius and aggregate to daily averages. The pre-downloaded NetCDF files (`boston_historical.nc` and `boston_ssp370.nc`) contain 3-hourly surface air temperature for the grid cell nearest to Boston. We aggregate 8 consecutive 3-hour periods to approximate daily averages, ignoring daylight saving time---a simplification typical in many bias correction applications. If you're interested in applying these methods to a different location, the `download_data.jl` script demonstrates how to extract CMIP6 data from Google Cloud Storage.

#block[
#callout(
body: 
[
Before starting the lab, uncomment the `Pkg.instantiate()` line in the first code block and run it to install all required packages. This will take a few minutes the first time. After installation completes, comment the line back out to avoid reinstalling on subsequent runs.

]
, 
title: 
[
Before Starting
]
, 
background_color: 
rgb("#f7dddc")
, 
icon_color: 
rgb("#CC1914")
, 
icon: 
fa-exclamation()
, 
body_background_color: 
white
)
]
= Lab Implementation
<lab-implementation>
== Package Setup
<package-setup>
#block[
```julia
using Pkg
lab_dir = dirname(@__FILE__)
Pkg.activate(lab_dir)
# Pkg.instantiate() # uncomment this the first time you run the lab to install packages, then comment it back
```

]
#block[
```julia
using CSV, CairoMakie, DataFrames, Dates, LaTeXStrings, NCDatasets, Statistics, StatsBase, TidierData
ENV["DATAFRAMES_ROWS"] = 6
CairoMakie.activate!()

# Constants for plotting
const MONTH_NAMES = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
const MONTH_LABELS = (1:12, MONTH_NAMES)
```

]
== Task 1: Load and Process Data
<task-1-load-and-process-data>
We begin by loading both observational and climate model data. The observational data from GHCN-Daily spans 1936-2024, while the CMIP6 data requires processing from 3-hourly to daily resolution. We'll use the overlapping historical period (1995-2014) to calibrate bias corrections.

#block[
#callout(
body: 
[
Load the Boston GHCN temperature data from `USW00014739.csv` and filter to years with at least 80% complete data. Then load the pre-downloaded CMIP6 NetCDF files and aggregate the 3-hourly temperature data to daily values for the historical (1995-2014) and near-term future (2020-2040) periods. The helper functions `load_cmip6_data()` and `aggregate_to_daily()` are provided below. Finally, visualize the annual cycle comparing observations vs the historical GCM simulation.

]
, 
title: 
[
Instructions
]
, 
background_color: 
rgb("#f7dddc")
, 
icon_color: 
rgb("#CC1914")
, 
icon: 
fa-exclamation()
, 
body_background_color: 
white
)
]
#block[
```julia
# Load and clean the observational data
data_path = joinpath(lab_dir, "USW00014739.csv")
df = @chain begin
    CSV.read(data_path, DataFrame)
    @mutate(
        TAVG = ifelse.(ismissing.(TMIN) .| ismissing.(TMAX), missing, (TMIN .+ TMAX) ./ 2),
    )
    @mutate(TAVG = TAVG / 10.0) # Convert to degrees C
    @rename(date = DATE)
    @mutate(year = year(date), month = month(date))
    @select(date, year, month, TAVG)
end

# Filter to years with at least 80% complete data
yearly_counts = @chain df begin
    @group_by(year)
    @summarize(n_obs = sum(!ismissing(TAVG)), n_total = n())
    @mutate(frac_complete = n_obs / n_total)
end

good_years_df = @chain yearly_counts begin
    @filter(frac_complete >= 0.8)
end
good_years = good_years_df.year

df_clean = @chain df begin
    @filter(year in !!good_years)
    dropmissing(:TAVG)
end
```

]
#block[
```julia
# Helper functions for CMIP6 data
"""
Load CMIP6 temperature data from a local NetCDF file and convert times to DateTime.
"""
function load_cmip6_data(file_path::String)
    ds = NCDataset(file_path)
    tas_data = ds["tas"][:]
    time_cf = ds["time"][:]
    close(ds)

    time_data = [DateTime(
        Dates.year(t), Dates.month(t), Dates.day(t),
        Dates.hour(t), Dates.minute(t), Dates.second(t)
    ) for t in time_cf]

    return tas_data, time_data
end

"""
Aggregate 3-hourly temperature data to daily averages.
"""
function aggregate_to_daily(tas_3hr, time_3hr)
    n_3hr_per_day = 8
    n_days = div(length(tas_3hr), n_3hr_per_day)

    daily_temp = Vector{Float64}(undef, n_days)
    daily_dates = Vector{Date}(undef, n_days)

    for i in 1:n_days
        idx_start = (i - 1) * n_3hr_per_day + 1
        idx_end = i * n_3hr_per_day
        daily_vals = collect(skipmissing(tas_3hr[idx_start:idx_end]))
        daily_temp[i] = isempty(daily_vals) ? NaN : mean(daily_vals)
        daily_dates[i] = Date(time_3hr[idx_start])
    end

    daily_temp_c = daily_temp .- 273.15  # Convert K to C
    return daily_temp_c, daily_dates
end

# Load and process CMIP6 data
hist_file = joinpath(lab_dir, "boston_historical.nc")
ssp370_file = joinpath(lab_dir, "boston_ssp370.nc")

tas_hist_3hr, time_hist_3hr = load_cmip6_data(hist_file)
tas_ssp370_3hr, time_ssp370_3hr = load_cmip6_data(ssp370_file)

# Process historical period (1995-2014)
hist_start = DateTime(1995, 1, 1)
hist_end = DateTime(2014, 12, 31, 23, 59, 59)
hist_idx = (time_hist_3hr .>= hist_start) .& (time_hist_3hr .<= hist_end)
tas_hist_daily, dates_hist_daily = aggregate_to_daily(tas_hist_3hr[hist_idx], time_hist_3hr[hist_idx])

# Process near-term period (2020-2040)
near_start = DateTime(2020, 1, 1)
near_end = DateTime(2040, 12, 31, 23, 59, 59)
near_idx = (time_ssp370_3hr .>= near_start) .& (time_ssp370_3hr .<= near_end)
tas_ssp370_near_daily, dates_ssp370_near_daily = aggregate_to_daily(tas_ssp370_3hr[near_idx], time_ssp370_3hr[near_idx])

# Create DataFrames
df_gcm_hist = DataFrame(
    date=dates_hist_daily, temp=tas_hist_daily,
    year=year.(dates_hist_daily), month=month.(dates_hist_daily)
)

df_ssp370_near = DataFrame(
    date=dates_ssp370_near_daily, temp=tas_ssp370_near_daily,
    year=year.(dates_ssp370_near_daily), month=month.(dates_ssp370_near_daily)
)

df_obs_hist = @chain df_clean begin
    @filter(year >= 1995 && year <= 2014)
end

# Compute monthly climatologies
obs_monthly = @chain df_obs_hist begin
    @group_by(month)
    @summarize(mean_temp = mean(TAVG))
end

gcm_monthly = @chain df_gcm_hist begin
    @group_by(month)
    @summarize(mean_temp = mean(temp))
end
```

]
```julia
let
    fig = Figure()
    ax = Axis(fig[1, 1],
        xlabel="Month",
        ylabel=L"Temperature ($^\circ$C)",
        title="Annual Cycle: Observations vs GCM Historical (1995-2014)",
        xticks=MONTH_LABELS)
    lines!(ax, obs_monthly.month, obs_monthly.mean_temp, linewidth=2, color=:steelblue, label="Observations")
    lines!(ax, gcm_monthly.month, gcm_monthly.mean_temp, linewidth=2, color=:coral, label="GCM Historical")
    axislegend(ax, position=:lt)
    fig
end
```

#figure([
#box(image("index_files/figure-typst/fig-gcm-obs-comparison-output-1.png", height: 5in, width: 7in))
], caption: figure.caption(
position: bottom, 
[
Annual cycle comparison between GHCN observations and GFDL-ESM4 historical simulation for Boston (1995-2014). The GCM shows a warm bias across most months.
]), 
kind: "quarto-float-fig", 
supplement: "Figure", 
)
<fig-gcm-obs-comparison>


== Task 2: Implement Delta Method
<task-2-implement-delta-method>
The delta method corrects the mean bias while preserving the climate model's projected change signal. We calculate a monthly bias correction based on the historical period (1995-2014) and apply it to future projections.

#block[
#callout(
body: 
[
Implement the additive delta method for temperature bias correction. For each calendar month $m$, calculate the mean bias: $Delta_m = macron(T)_(upright("hist") \, m)^(upright("GCM")) - macron(T)_(upright("hist") \, m)^(upright("obs"))$. Then apply the correction to future values: $T_(upright("fut"))^(upright("corr")) (d \, m \, y) = T_(upright("fut"))^(upright("GCM")) (d \, m \, y) - Delta_m$.

Follow these steps:

+ Calculate the monthly mean bias by grouping both `df_gcm_hist` and `df_obs_hist` by month, computing their means, joining them, and computing the difference.
+ Create a bar plot visualizing the monthly bias.
+ Write a function `apply_delta_method(gcm_temps, gcm_dates, monthly_bias_df)` that applies the bias correction to a vector of temperatures and dates.
+ Apply your function to the near-term data and add a new column `temp_delta` to `df_ssp370_near`.
+ Create monthly climatologies for both raw and delta-corrected near-term data.
+ Visualize the annual cycle showing historical observations, raw GCM near-term, and delta-corrected near-term temperatures.

]
, 
title: 
[
Instructions
]
, 
background_color: 
rgb("#f7dddc")
, 
icon_color: 
rgb("#CC1914")
, 
icon: 
fa-exclamation()
, 
body_background_color: 
white
)
]
#block[
```julia
# Step 1 - Calculate the mean difference: Δ_m = mean(GCM_hist) - mean(Obs_hist)
gcm_monthly = @chain df_gcm_hist begin
    @group_by(month)
    @summarize(gcm_mean = mean(temp))
    @arrange(month)
end

obs_monthly = @chain df_obs_hist begin
    @group_by(month)
    @summarize(obs_mean = mean(TAVG))
    @arrange(month)
end

monthly_bias = leftjoin(gcm_monthly, obs_monthly, on = :month)
monthly_bias.bias = monthly_bias.gcm_mean .- monthly_bias.obs_mean
monthly_bias = sort(monthly_bias, :month)

@assert all(sort(unique(monthly_bias.month)) .== collect(1:12)) "Monthly bias missing some months"

# Print the GEM and Obs data as dataframe
for row in eachrow(monthly_bias)
    println("Month: $(row.month), GCM mean = $(round(row.gcm_mean, digits=2)), ",
            "Obs mean = $(round(row.obs_mean, digits=2)), ",
            "Bias (GCM - Obs) = $(round(row.bias, digits=2))")
end

```

#block[
```
Month: 1, GCM mean = -3.18, Obs mean = -1.09, Bias (GCM - Obs) = -2.09
Month: 2, GCM mean = -1.7, Obs mean = 0.02, Bias (GCM - Obs) = -1.73
Month: 3, GCM mean = 2.01, Obs mean = 3.79, Bias (GCM - Obs) = -1.78
Month: 4, GCM mean = 6.41, Obs mean = 9.34, Bias (GCM - Obs) = -2.93
Month: 5, GCM mean = 12.44, Obs mean = 14.59, Bias (GCM - Obs) = -2.14
Month: 6, GCM mean = 17.86, Obs mean = 19.94, Bias (GCM - Obs) = -2.08
Month: 7, GCM mean = 21.14, Obs mean = 23.34, Bias (GCM - Obs) = -2.2
Month: 8, GCM mean = 20.55, Obs mean = 22.54, Bias (GCM - Obs) = -1.98
Month: 9, GCM mean = 17.18, Obs mean = 18.69, Bias (GCM - Obs) = -1.51
Month: 10, GCM mean = 11.56, Obs mean = 12.69, Bias (GCM - Obs) = -1.14
Month: 11, GCM mean = 5.15, Obs mean = 7.06, Bias (GCM - Obs) = -1.92
Month: 12, GCM mean = 0.52, Obs mean = 2.05, Bias (GCM - Obs) = -1.53
```

]
]
The GCM model underestimate so the bias is negative. Need to adjust the model.

```julia
# Step 2 - Barplot for monthly bias（GCM - Obs，1995–2014）
fig = Figure()
ax  = Axis(fig[1, 1];
    title = "Monthly Mean Bias (GCM − Obs), 1995–2014",
    xlabel = "Month",
    ylabel = L"Temperature Bias ($^\circ$C)",
    xticks = MONTH_LABELS
)

barplot!(ax, monthly_bias.month, monthly_bias.bias, strokewidth = 0)
hlines!(ax, 0; color = :black, linestyle = :dash, linewidth = 1)

fig
```

#box(image("index_files/figure-typst/cell-8-output-1.png", height: 5in, width: 7in))

```julia
# Stpe 3 - Apply the delta-method bias correction: T_corr = T_gcm - Δ_m.


# gcm_temps::AbstractVector: Future-period GCM temperatures
# gcm_dates::AbstractVector{<:Date}: Dates corresponding to each temperature
# monthly_bias_df::DataFrame: Must contain columns `month` and `bias`

function apply_delta_method(gcm_temps, gcm_dates, monthly_bias_df)
    corrected_temps = similar(gcm_temps)  # preserve element type & length

    # Create a dictionary for fast month-to-bias lookup
    bias_map = Dict(monthly_bias_df.month[i] => monthly_bias_df.bias[i] 
                    for i in 1:nrow(monthly_bias_df))

    # Apply correction
    for i in eachindex(gcm_temps)
        m = month(gcm_dates[i])
        Δm = bias_map[m]  # lookup bias for that month
        corrected_temps[i] = gcm_temps[i] - Δm
    end

    return corrected_temps
end
```

```
apply_delta_method (generic function with 1 method)
```

```julia
# Step 4 - Apply the correction

df_ssp370_near.temp_delta = apply_delta_method(
    df_ssp370_near.temp,
    df_ssp370_near.date,
    monthly_bias
)
first(df_ssp370_near, 5)
```

#table(
  columns: 6,
  align: (right,left,right,right,right,right,),
  table.header(table.cell(align: right)[#set text(weight: "bold"); Row], table.cell(align: left)[date], table.cell(align: left)[temp], table.cell(align: left)[year], table.cell(align: left)[month], table.cell(align: left)[temp\_delta],
    table.cell(align: right)[#set text(weight: "bold"); ], table.cell(align: left)[Date], table.cell(align: left)[Float64], table.cell(align: left)[Int64], table.cell(align: left)[Int64], table.cell(align: left)[Float64],),
  table.hline(),
  table.cell(align: right)[#set text(weight: "bold"); 1], table.cell(align: left)[2020-01-01], table.cell(align: right)[-7.3856], table.cell(align: right)[2020], table.cell(align: right)[1], table.cell(align: right)[-5.2931],
  table.cell(align: right)[#set text(weight: "bold"); 2], table.cell(align: left)[2020-01-02], table.cell(align: right)[-7.71573], table.cell(align: right)[2020], table.cell(align: right)[1], table.cell(align: right)[-5.62324],
  table.cell(align: right)[#set text(weight: "bold"); 3], table.cell(align: left)[2020-01-03], table.cell(align: right)[-3.07117], table.cell(align: right)[2020], table.cell(align: right)[1], table.cell(align: right)[-0.978678],
  table.cell(align: right)[#set text(weight: "bold"); 4], table.cell(align: left)[2020-01-04], table.cell(align: right)[0.409509], table.cell(align: right)[2020], table.cell(align: right)[1], table.cell(align: right)[2.502],
  table.cell(align: right)[#set text(weight: "bold"); 5], table.cell(align: left)[2020-01-05], table.cell(align: right)[3.31857], table.cell(align: right)[2020], table.cell(align: right)[1], table.cell(align: right)[5.41106],
)
#block[
```julia
# Step 5 - Compute monthly climatologies
near_monthly_raw = @chain df_ssp370_near begin
    @group_by(month)
    @summarize(mean_temp = mean(temp))
    @arrange(month)
end

near_monthly_delta = @chain df_ssp370_near begin
    @group_by(month)
    @summarize(mean_temp = mean(temp_delta))
    @arrange(month)
end

@assert all(near_monthly_raw.month .== collect(1:12))
@assert all(near_monthly_delta.month .== collect(1:12))

println(names(near_monthly_raw))
println(near_monthly_raw)

println(names(near_monthly_delta))
println(near_monthly_delta)
```

#block[
```
["month", "mean_temp"]
12×2 DataFrame
 Row │ month  mean_temp  
     │ Int64  Float64    
─────┼───────────────────
   1 │     1  -1.98899
   2 │     2  -0.0236623
   3 │     3   2.83692
   4 │     4   7.77089
   5 │     5  12.8063
   6 │     6  17.9288
   7 │     7  21.6318
   8 │     8  21.622
   9 │     9  17.9171
  10 │    10  12.3487
  11 │    11   5.9952
  12 │    12   0.985578
["month", "mean_temp"]
12×2 DataFrame
 Row │ month  mean_temp 
     │ Int64  Float64   
─────┼──────────────────
   1 │     1   0.103505
   2 │     2   1.70296
   3 │     3   4.62097
   4 │     4  10.7028
   5 │     5  14.9487
   6 │     6  20.0121
   7 │     7  23.8329
   8 │     8  23.6061
   9 │     9  19.4281
  10 │    10  13.4868
  11 │    11   7.91067
  12 │    12   2.51416
```

]
]
```julia
# Step 6 - Creat comparison plot
rename!(obs_monthly, :obs_mean => :mean_temp)

fig = Figure()
ax  = Axis(fig[1, 1];
    title  = "Delta Method: Annual Cycle for SSP3-7.0 Near-Term (2020–2040)",
    xlabel = "Month",
    ylabel = L"Temperature ($^\circ$C)",
    xticks = MONTH_LABELS
)

# Historical Obs
lines!(ax, obs_monthly.month, obs_monthly.mean_temp;
    linewidth = 2, color = :steelblue, label = "Historical Obs")

# Near monthly raw
lines!(ax, near_monthly_raw.month, near_monthly_raw.mean_temp;
    linewidth = 2, color = :coral, label = "GCM Raw")

# Near monthly delta
lines!(ax, near_monthly_delta.month, near_monthly_delta.mean_temp;
    linewidth = 2, color = :green, linestyle = :dash, label = "Delta Corrected")

axislegend(ax, position = :lt)
fig
```

#box(image("index_files/figure-typst/cell-12-output-1.png", height: 5in, width: 7in))

== Task 3: Implement Quantile-Quantile Mapping
<task-3-implement-quantile-quantile-mapping>
Unlike the delta method, QQ-mapping transforms the entire probability distribution of model output to match observations. For each value in the future model output, we find its percentile in the historical model distribution, then map it to the same percentile in the historical observed distribution.

#block[
#callout(
body: 
[
Implement QQ-mapping to correct the full distribution of temperature values.

Follow these steps:

+ Extract the observed and GCM historical temperatures as vectors (`obs_hist_temps` and `gcm_hist_temps`).
+ Use `ecdf()` from StatsBase.jl to fit empirical cumulative distribution functions to both datasets.
+ Write a function `apply_qqmap_empirical(gcm_temps, ecdf_gcm, obs_hist_temps)` that:
  - For each temperature value, finds its percentile in the GCM CDF
  - Clamps the percentile to the range \[0.001, 0.999\] to avoid extrapolation
  - Maps to the same percentile in the observed distribution using `quantile()`
+ Apply your function to the near-term data and add a new column `temp_qqmap` to `df_ssp370_near`.
+ Create a QQ-plot comparing observed and GCM quantiles for the historical period, showing the 1:1 line.

#block[
#callout(
body: 
[
- Use `ecdf_gcm(value)` to get the percentile of a value in the GCM distribution
- Use `quantile(obs_hist_temps, p)` to get the value at percentile `p` in the observed distribution
- The `clamp(x, low, high)` function constrains `x` to the range \[low, high\]

]
, 
title: 
[
Hints
]
, 
background_color: 
rgb("#ccf1e3")
, 
icon_color: 
rgb("#00A047")
, 
icon: 
fa-lightbulb()
, 
body_background_color: 
white
)
]
]
, 
title: 
[
Instructions
]
, 
background_color: 
rgb("#f7dddc")
, 
icon_color: 
rgb("#CC1914")
, 
icon: 
fa-exclamation()
, 
body_background_color: 
white
)
]
#block[
```julia
# Step 1 - Extract temperature vectors as Float64 (historical period)

# Obs (1995 - 2014)：df_obs_hist.TAVG -> Float64
obs_hist_temps = Float64.(collect(skipmissing(df_obs_hist.TAVG))) # skip missing

# GCM（历史期 1995–2014）：df_gcm_hist.temp -> Float64 
gcm_hist_temps = Float64.(collect(skipmissing(df_gcm_hist.temp))) # skip missing


obs_hist_temps = filter(isfinite, obs_hist_temps)
gcm_hist_temps = filter(isfinite, gcm_hist_temps)


@assert length(obs_hist_temps) > 0 "obs_hist_temps is empty"
@assert length(gcm_hist_temps) > 0 "gcm_hist_temps is empty"
println("Obs hist n=", length(obs_hist_temps), " | GCM hist n=", length(gcm_hist_temps))
```

#block[
```
Obs hist n=7305 | GCM hist n=7300
```

]
]
#block[
```julia
# Step 2 - Fit empirical CDFs

ecdf_obs = ecdf(obs_hist_temps)
ecdf_gcm = ecdf(gcm_hist_temps)

println("ECDF check — obs at 0°C:   ", ecdf_obs(0.0))
println("ECDF check — gcm at 0°C:   ", ecdf_gcm(0.0))
```

#block[
```
ECDF check — obs at 0°C:   0.13853524982888432
ECDF check — gcm at 0°C:   0.1917808219178082
```

]
]
```julia
# Step 3 - Write the mapping function

"""
Apply empirical Quantile-Quantile (QQ) mapping correction.

Arguments:
- gcm_temps::Vector{Float64}: Future-period GCM temperature data to correct
- ecdf_gcm: Empirical CDF function of historical GCM temps (from StatsBase.ecdf)
- obs_hist_temps::Vector{Float64}: Historical observed temps for quantile mapping

Returns:
- corrected_temps::Vector{Float64}: QQ-mapped corrected temperatures
"""
function apply_qqmap_empirical(gcm_temps, ecdf_gcm, obs_hist_temps)
    corrected_temps = similar(gcm_temps)  # same type and size as input

    for i in eachindex(gcm_temps)
        # Step 1: Get percentile from historical GCM distribution
        p = ecdf_gcm(gcm_temps[i])

        # Step 2: Clamp percentile to avoid extreme tail errors
        p_clamped = clamp(p, 0.001, 0.999)

        # Step 3: Map percentile into observed historical distribution
        corrected_temps[i] = quantile(obs_hist_temps, p_clamped)
    end

    return corrected_temps
end

```

```
Main.Notebook.apply_qqmap_empirical
```

```julia
# Step 4 - Apply the correction

df_ssp370_near.temp_qqmap = apply_qqmap_empirical(
    Float64.(df_ssp370_near.temp),  # Future GCM temperatures
    ecdf_gcm,                       # Historical GCM ECDF
    obs_hist_temps                 # Historical OBS distribution
)
first(df_ssp370_near, 5)
```

#table(
  columns: 7,
  align: (right,left,right,right,right,right,right,),
  table.header(table.cell(align: right)[#set text(weight: "bold"); Row], table.cell(align: left)[date], table.cell(align: left)[temp], table.cell(align: left)[year], table.cell(align: left)[month], table.cell(align: left)[temp\_delta], table.cell(align: left)[temp\_qqmap],
    table.cell(align: right)[#set text(weight: "bold"); ], table.cell(align: left)[Date], table.cell(align: left)[Float64], table.cell(align: left)[Int64], table.cell(align: left)[Int64], table.cell(align: left)[Float64], table.cell(align: left)[Float64],),
  table.hline(),
  table.cell(align: right)[#set text(weight: "bold"); 1], table.cell(align: left)[2020-01-01], table.cell(align: right)[-7.3856], table.cell(align: right)[2020], table.cell(align: right)[1], table.cell(align: right)[-5.2931], table.cell(align: right)[-5.75],
  table.cell(align: right)[#set text(weight: "bold"); 2], table.cell(align: left)[2020-01-02], table.cell(align: right)[-7.71573], table.cell(align: right)[2020], table.cell(align: right)[1], table.cell(align: right)[-5.62324], table.cell(align: right)[-6.1],
  table.cell(align: right)[#set text(weight: "bold"); 3], table.cell(align: left)[2020-01-03], table.cell(align: right)[-3.07117], table.cell(align: right)[2020], table.cell(align: right)[1], table.cell(align: right)[-0.978678], table.cell(align: right)[-1.95],
  table.cell(align: right)[#set text(weight: "bold"); 4], table.cell(align: left)[2020-01-04], table.cell(align: right)[0.409509], table.cell(align: right)[2020], table.cell(align: right)[1], table.cell(align: right)[2.502], table.cell(align: right)[2.5],
  table.cell(align: right)[#set text(weight: "bold"); 5], table.cell(align: left)[2020-01-05], table.cell(align: right)[3.31857], table.cell(align: right)[2020], table.cell(align: right)[1], table.cell(align: right)[5.41106], table.cell(align: right)[5.8],
)
```julia
# Step 5 - Creat the QQ plot

percentiles = collect(0.01:0.01:0.99) # Define precentiles

# Compute two quantiles
obs_quantiles      = [quantile(obs_hist_temps,  p) for p in percentiles]
gcm_hist_quantiles = [quantile(gcm_hist_temps, p) for p in percentiles]

# Scatter plot 
fig = Figure()
ax  = Axis(fig[1, 1];
    title  = "QQ Plot: GCM Historical vs Observations (1995–2014)",
    xlabel = L"Observed Quantiles ($^\circ$C)",
    ylabel = L"GCM Quantiles ($^\circ$C)",
    aspect = DataAspect() # plot square
)

scatter!(ax, obs_quantiles, gcm_hist_quantiles; markersize=6)
lines!(ax, [-20, 40], [-20, 40]; linestyle=:dash, color=:black)  # 1:1 reference line

fig
```

#box(image("index_files/figure-typst/cell-17-output-1.png", height: 5in, width: 7in))

== Task 4: Compare Methods
<task-4-compare-methods>
Now we compare the delta method and QQ-mapping approaches to understand their strengths and limitations.

#block[
#callout(
body: 
[
Create a single figure showing monthly mean temperature for the near-term period (2020-2040) using all four approaches:

+ Compute the monthly mean for QQ-mapped temperatures (`near_monthly_qqmap`)
+ Create a figure with 4 lines:
  - Historical observations (`obs_monthly`)
  - Raw GCM near-term (`near_monthly_raw`)
  - Delta-corrected near-term (`near_monthly_delta`)
  - QQ-mapped near-term (`near_monthly_qqmap`)

After creating the plot, examine it carefully and consider:

- How do the two correction methods differ?
- What does the QQ-mapped line reveal about this method's treatment of the warming signal?

]
, 
title: 
[
Instructions
]
, 
background_color: 
rgb("#f7dddc")
, 
icon_color: 
rgb("#CC1914")
, 
icon: 
fa-exclamation()
, 
body_background_color: 
white
)
]
```julia
# Step 1 - Monthy mean for QQ-mapped temepratures

near_monthly_qqmap = @chain df_ssp370_near begin
    @group_by(month)
    @summarize(mean_temp = mean(temp_qqmap))
    @arrange(month)
end
```

#table(
  columns: 3,
  align: (right,right,right,),
  table.header(table.cell(align: right)[#set text(weight: "bold"); Row], table.cell(align: left)[month], table.cell(align: left)[mean\_temp],
    table.cell(align: right)[#set text(weight: "bold"); ], table.cell(align: left)[Int64], table.cell(align: left)[Float64],),
  table.hline(),
  table.cell(align: right)[#set text(weight: "bold"); 1], table.cell(align: right)[1], table.cell(align: right)[-0.374928],
  table.cell(align: right)[#set text(weight: "bold"); 2], table.cell(align: right)[2], table.cell(align: right)[1.75688],
  table.cell(align: right)[#set text(weight: "bold"); 3], table.cell(align: right)[3], table.cell(align: right)[4.83506],
  table.cell(align: right)[⋮], table.cell(align: right)[⋮], table.cell(align: right)[⋮],
  table.cell(align: right)[#set text(weight: "bold"); 10], table.cell(align: right)[10], table.cell(align: right)[14.1417],
  table.cell(align: right)[#set text(weight: "bold"); 11], table.cell(align: right)[11], table.cell(align: right)[8.04862],
  table.cell(align: right)[#set text(weight: "bold"); 12], table.cell(align: right)[12], table.cell(align: right)[2.77888],
)
```julia
# Step 2 - Create the formatted comparison figure

fig = Figure()
ax = Axis(
    fig[1, 1];
    title = "Annual Cycle Comparison: All Methods (Near-term 2020–2040)",
    xlabel = "Month",
    ylabel = L"Temperature ($^\circ$C)",
    xticks = MONTH_LABELS
)

# Historical Observations (black, solid, width = 2.5)
lines!(ax,
    obs_monthly.month, obs_monthly.mean_temp;
    color = :black, linewidth = 2.5, linestyle = :solid, label = "Historical Obs"
)

# Raw GCM near-term (coral/red, solid)
lines!(ax,
    near_monthly_raw.month, near_monthly_raw.mean_temp;
    color = :coral, linewidth = 2, linestyle = :solid, label = "Raw GCM"
)

# Delta corrected (green, dashed)
lines!(ax,
    near_monthly_delta.month, near_monthly_delta.mean_temp;
    color = :green, linewidth = 2, linestyle = :dash, label = "Delta"
)

# QQ-mapped (purple, dotted)
lines!(ax,
    near_monthly_qqmap.month, near_monthly_qqmap.mean_temp;
    color = :purple, linewidth = 2, linestyle = :dot, label = "QQ-map"
)

axislegend(ax, position = :lt)
fig

```

#box(image("index_files/figure-typst/cell-19-output-1.png", height: 5in, width: 7in))

The delta method applies a fixed mean adjustment to each month which results in a uniform shift of the future temperature cycle while preserving the model projected warming pattern. In contrast, the QQ mapping method remaps each future value according to its empirical percentile, which causes a redistribution of the full temperature distribution. The redistribution reduces the magnitude of projected warming because it constrains future values to align more closely with the historical observed distribution. The elevated summer peak in the QQ mapped line indicates that this approach moderates extremes instead of allowing the stronger warming predicted by the raw model.

The QQ mapped curve lies closer to the historical observations in several months, like spring time when the tempwerature is soft, which implies stronger adherence to the reference climatology. The tighter alignment can improve representation of present day conditions but it may artificially weaken future warming signals. Such dampening of the climate change signal can underestimate projected risks in impact assessments and may lead to insufficient adaptation or resilience planning when future extremes are critical to decision making.

== Task 5: Reflection
<task-5-reflection>
Finally, we reflect on the assumptions, limitations, and appropriate use cases for each bias correction method.

#block[
#callout(
body: 
[
Write brief responses (2-3 sentences each) to the following questions:

+ #strong[Method selection:] If you were providing climate data to support a decision about urban heat management in Boston, which bias correction method would you recommend and why?

+ #strong[Appropriateness of QQ-mapping:] In the #cite(<ines_biascorrection:2006>, form: "prose") paper we discussed, QQ-mapping was used to correct both the #emph[frequency] and #emph[intensity] of rainfall for crop modeling. For temperature data, does it make sense to correct the full distribution? Why or why not? When might the delta method be more appropriate than QQ-mapping?

]
, 
title: 
[
Instructions
]
, 
background_color: 
rgb("#f7dddc")
, 
icon_color: 
rgb("#CC1914")
, 
icon: 
fa-exclamation()
, 
body_background_color: 
white
)
]
=== Method Selection for Urban Heat Management
<method-selection-for-urban-heat-management>
#emph[For urban heat management, I would recommend QQ-mapping because it corrects the full distribution of temperatures rather than only adjusting monthly means. The QQ-plot demonstrates that QQ-mapping better aligns the model with the observed distribution across all percentiles. The method could provides a more reliable basis for evaluating thresholds relevant to public health and urban cooling strategies.]

=== Appropriateness of QQ-Mapping for Temperature
<appropriateness-of-qq-mapping-for-temperature>
#emph[In the study by Ines and Hansen (2006), QQ-mapping was essential for rainfall because impact models such as crop simulations are sensitive to both event frequency and intensity, especially under highly skewed precipitation distributions. For temperature, the distribution is typically smoother and less intermittent. Full distribution mapping may sometimes overconstrain future projections and dampen climate change signals that policymakers must consider. In cases where preserving the magnitude of projected warming is a priority, the delta method may be more appropriate because it maintains the climate model's change signal while still correcting systematic mean bias. ]

= References
<references>
#block[
#block[
```julia

open("references.bib", "w") do f
    write(f, """
@article{_biascorrection:2006,
  author  = {Ines, A. V. M. and Hansen, J. W.},
  title   = {Bias Correction of Daily GCM Rainfall for Crop Simulation Studies},
  journal = {Agricultural and Forest Meteorology},
  year    = {2006},
  volume  = {138},
  number  = {1},
  pages   = {44--53},
  doi     = {10.1016/j.agrformet.2006.03.009}
}
""")
end

println("Bib file exist")
```

#block[
```
Bib file exist
```

]
]
] <refs>




#bibliography("references.bib")

