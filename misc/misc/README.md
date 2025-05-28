# Get Started

We provide a script (tested under python 3.10.9) to quickly analyze (usually a few seconds) the output of centroAnno. This script only needs the package matplotlib.
```bash
# Install
pip install matplotlib

# You can quickly analyze (usually a few seconds) the output of centroAnno
python cautils.py $centroAnno_output_dir $your_analysis_dir

# Example command (The output `tmp1` should be the same as the `tmp` we provided)
python cautils.py test tmp1
```

# Output

`repeat_regions.bed`: Tandem repeat regions greater than a certain length (1000bp). Column 1: Sequence Name, Column 2: Start Position, Column 3: End Position, Column 4: Most Frequent Tandem Repeat Length.

`top10_repeats.bed`: Top 10 tandem repeat unit lengths. Column 1: Sequence Name, Column 2: Unit Length, Column 3: Span Length.

`HORs.bed`: High confidence higher-order tandem repeats (HORs). Column 1: Sequence Name, Column 2: HOR Name,  Column 3: Start Position, Column 4: End Position, Column 5: The Number of Monomers Contained In HOR, Column 6: HOR Length, Column 7: Span Length.

`mono.svg`: Visualization of tandem repeat units.

`hor.svg`: Visualization of HORs.
