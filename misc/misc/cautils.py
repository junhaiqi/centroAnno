# A script for centroAnno's output analysis

import sys
import matplotlib.pyplot as plt
import os
from collections import Counter
import random
# matplotlib.rcParams['font.family'] = 'Times New Roman'
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from matplotlib.colorbar import ColorbarBase
# from matplotlib import colormaps
# cmap = colormaps.get_cmap("viridis")
cmap = plt.colormaps["viridis"]
from mpl_toolkits.axes_grid1 import make_axes_locatable

repeat_region_minlen = 1000
hor_confidence_cutoff = 0.9
fig_fontsize = 18
plt.rcParams['font.size'] = fig_fontsize

def find_mode(lst):
    if not lst:
        return None, 0 
    counter = Counter(lst)  
    mode, count = counter.most_common(1)[0]  
    return mode, count

def sort_bed_file(input_file, output_file):
    cmd = f'./bedtools sort -i {input_file} >{output_file}'
    os.system( cmd )
    
def draw_mono_fig(bed_file, out_fig):
    data = []
    with open(bed_file) as f:
        lines = f.readlines()
        for line in lines:
            info = line.strip('\n').split('\t')
            this_data = [info[0], int(info[1]), int( info[2] ), int( info[3] )]
            data.append( this_data )
            
    hor_lengths = [row[3] for row in data]
    min_len, max_len = min(hor_lengths), max(hor_lengths)
    norm = mcolors.Normalize(vmin=min_len, vmax=max_len)
    if min_len == max_len:
        norm = mcolors.Normalize(vmin=min_len, vmax=max_len + 1)
    all_seqs = sorted(set(row[0] for row in data))
    seq_y_offset = {seq: i * 20 for i, seq in enumerate(all_seqs)}
    fig, ax = plt.subplots(figsize=(24, 9))
    for row in data:
        seq, start, end, hor_len = row
        ypos = seq_y_offset[seq]
        color = cmap(norm(hor_len))
        ax.broken_barh([(start, end - start)], (ypos, 8), facecolors=color)
        
    for seq in all_seqs:
        ypos = seq_y_offset[seq]
        ax.text(0, ypos + 4, seq, va='center', ha='right', fontsize=12, fontweight='bold')
        
    ax.set_ylim(-10, max(seq_y_offset.values()) + 20)
    ax.set_xlim(0, max(row[2] for row in data) + 300)
    ax.set_xlabel("Genomic Coordinate")
    ax.set_yticks([])
    ax.set_title("Tandem Repeat Unit Visualization")
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.1)
    cb = ColorbarBase(cax, cmap=cmap, norm=norm, orientation='vertical')
    cb.set_label("Tandem Repeat Unit Length (bp)")
    plt.tight_layout()
    plt.savefig(out_fig, bbox_inches='tight')
    
    
def draw_hor_fig(bed_file, out_fig):
    data = []
    with open(bed_file) as f:
        lines = f.readlines()
        for line in lines:
            info = line.strip('\n').split('\t')
            # this_data = [info[0], info[1], int( info[2] ), int( info[3] ), int( info[4] ), int( info[5] ), int( info[6] )]
            this_data = [info[0], info[1], int( info[2] ), int( info[3] ), int( info[5] ) ]
            data.append( this_data )

    hor_lengths = [row[4] for row in data]
    min_len, max_len = min(hor_lengths), max(hor_lengths)

    norm = mcolors.Normalize(vmin=min_len, vmax=max_len)
    if min_len == max_len:
        norm = mcolors.Normalize(vmin=min_len, vmax=max_len + 1)
    all_seqs = sorted(set(row[0] for row in data))
    seq_y_offset = {seq: i * 20 for i, seq in enumerate(all_seqs)}

    fig, ax = plt.subplots(figsize=(24, 9))
    
    for row in data:
        seq, hor_type, start, end, hor_len = row
        ypos = seq_y_offset[seq]
        color = cmap(norm(hor_len))
        ax.broken_barh([(start, end - start)], (ypos, 8), facecolors=color)

    
    for seq in all_seqs:
        ypos = seq_y_offset[seq]
        ax.text(0, ypos + 4, seq, va='center', ha='right', fontsize=12, fontweight='bold')
        
    ax.set_ylim(-10, max(seq_y_offset.values()) + 20)
    ax.set_xlim(0, max(row[3] for row in data) + 300)
    ax.set_xlabel("Genomic Coordinate")
    ax.set_yticks([])
    ax.set_title("HOR Visualization")
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.1)
    cb = ColorbarBase(cax, cmap=cmap, norm=norm, orientation='vertical')
    cb.set_label("HOR Unit Length (bp)")
    plt.tight_layout()
    plt.savefig(out_fig, bbox_inches='tight')

def main(centroAnno_output_dir, analysis_out_dir):
    if os.path.exists(centroAnno_output_dir) == False:
        print(f"{centroAnno_output_dir} does not exist!")
        exit(-1)
    
    if os.path.exists(analysis_out_dir) == False:
        os.makedirs(analysis_out_dir)
    
    mono_dem_file = ""
    hor_dem_file = ""
    
    out_put_mono_bed = f'{analysis_out_dir}/repeat_regions.bed'
    out_put_mono_sorted_bed = f'{analysis_out_dir}/repeat_regions_sorted.bed'
    out_top10_monos_file = f'{analysis_out_dir}/top10_repeats.bed'
    out_high_confidence_hor_file = f'{analysis_out_dir}/HORs.bed'
    out_mono_fig = f'{analysis_out_dir}/mono.svg'
    out_hor_fig = f'{analysis_out_dir}/hor.svg'
    
    for item in os.listdir( centroAnno_output_dir ):
        if "_decomposedResult.csv" in item:
            mono_dem_file = f'{centroAnno_output_dir}/{item}'
        elif "_horDecomposedResult.csv" in item:
            hor_dem_file = f'{centroAnno_output_dir}/{item}'
            
    ########################## get repeat regions ##########################
    w_file = open(out_put_mono_bed, 'w')
    with open( mono_dem_file ) as f:
        lines = f.readlines()
        region_st = int(lines[1].strip('\n').split(',')[2])
        region_ed = int(lines[1].strip('\n').split(',')[3])
        rep_len_list = []
        for line in lines:
            if 'name' in line:
                continue
            info = line.strip('\n').split(',')
            st = int(info[2])
            ed = int(info[3])
            # iden = float(info[-2])
            if abs(region_ed - st) < 100:
                region_ed = ed
                rep_len_list.append( int( info[-1] ) )
            else:
                rep_len = find_mode( rep_len_list )[0]
                if abs(region_ed - region_st) >= repeat_region_minlen:
                    w_file.write(f'{info[0]}\t{region_st}\t{region_ed}\t{rep_len}\n')
                rep_len_list = []
                region_st = st
                region_ed = ed
        rep_len = find_mode( rep_len_list )[0]
        if abs(region_ed - region_st) >= repeat_region_minlen:
            w_file.write(f'{info[0]}\t{region_st}\t{region_ed}\t{rep_len}\n')
    w_file.close()
    sort_bed_file(out_put_mono_bed, out_put_mono_sorted_bed)
    ########################## get repeat regions ##########################
    
    ########################## get top10 repeats ##########################
    repeat_dict = {}
    with open( mono_dem_file ) as f:
        lines = f.readlines()
        for line in lines:
            if 'name' in line:
                continue
            info = line.strip('\n').split(',')
            rep_len = int( info[-1] )
            name = info[0]
            if name not in repeat_dict:
                repeat_dict[name] = {}
            key = str(rep_len) + 'bp'
            rep_base_num = int( info[3] ) - int( info[2] )
            if key in repeat_dict[name]:
                repeat_dict[name][key] += rep_base_num
            else:
                repeat_dict[name][key] = rep_base_num
                
    w_file = open( out_top10_monos_file, 'w')
    for name in repeat_dict:
        sub_repeat_dict = sorted(repeat_dict[name].items(), key=lambda x: x[1], reverse=True)
        top_num = min( 10, len( sub_repeat_dict ) )
        for i in range(top_num):
            w_file.write(f'{name}\t{sub_repeat_dict[i][0]}\t{sub_repeat_dict[i][1]}\n')
    w_file.close()
    ########################## get top10 repeats ##########################
    
    ########################## get HORs ##########################
    hor_dict = {}
    with open( hor_dem_file ) as f:
        lines = f.readlines()
        for line in lines:
            if 'name' in line:
                continue
            info = line.strip('\n').split(',')
            name = info[0]
            if name not in hor_dict:
                hor_dict[name] = {}
            hor_name = info[1]
            iden = float( info[-3] )
            st = int( info[2] )
            ed = int( info[3] )
            hor_len = int( info[-2] )
            hor_mon_num = hor_name.count('_') + 1
            if iden >= hor_confidence_cutoff and hor_mon_num > 1:
                if hor_name in hor_dict[name]:
                    if abs(hor_dict[name][ hor_name ][-1][1] - st) < 100:
                        hor_dict[name][ hor_name ][-1][1] = ed
                    else:
                        hor_dict[name][ hor_name ].append( [st, ed, hor_mon_num, hor_len] )
                else:
                    hor_dict[name][ hor_name ] = [ [st, ed, hor_mon_num, hor_len] ]
                                    
    w_file = open( out_high_confidence_hor_file, 'w')
    for key in hor_dict:
        for key1 in hor_dict[key]:
            mono_type_num = key1.count('_') + 1
            if mono_type_num > 1:
                for region in hor_dict[key][key1]:
                    span_len = region[1] - region[0]
                    if span_len >= region[3]:
                        w_file.write(f'{key}\t{key1}\t{ region[0] }\t{ region[1] }\t{ region[2]}\t{ region[3] }\t{span_len}\n')
    w_file.close()
    ########################## get HORs ##########################
    
    ########################## Visualization ##########################
    draw_mono_fig(out_put_mono_sorted_bed, out_mono_fig)
    draw_hor_fig(out_high_confidence_hor_file, out_hor_fig)
    ########################## Visualization ##########################

if __name__ == "__main__":
    if len( sys.argv ) != 3:
        print("*** A script for centroAnno's output analysis ***")
        print(f"Usage: python {sys.argv[0]} $centroAnno_output_dir $your_analysis_dir")
    else:
        main(sys.argv[1], sys.argv[2])
