#Need to create dictionary from input wig files (already normalized rpm from STAR)
#Then simply divide dictionaries from different samples and output to a new wig file output
#Also output a tab-delimeted file for MA plots with Chromosome\tposition\tavgreads\tfoldchange

import sys

#Read in two wig files from user and parse into nested dictionary containing 25S, 18S, 5.8S, 5S rRNA for each sample
#my_dict[expt1]["18S"][30] --> chooses expt 1 dictionary, 18S dictionary, then item 30
#mydict_expt1 = {'18S': {'1':0.13, '2':0.25}, '25S': {...}}

def parse_wig(filename):

    mydict = {}
    for line in open(filename):
        if "variableStep" in line:
            chrom = line.strip().split('=')[-1]
            mydict[chrom] = {}

        else:
            ll = line.strip().split('\t')
            mydict[chrom][int(ll[0])] = float(ll[1])

    return mydict


def make_ratio(numerator_dict, denominator_dict):

    ratio_dict = {}
    for chrom in numerator_dict:

        ratio_dict[chrom] ={}
        for pos in numerator_dict[chrom]:
            if pos in denominator_dict[chrom]:
                ratio_dict[chrom][pos] = numerator_dict[chrom][pos] / denominator_dict[chrom][pos]
            else:
                ratio_dict[chrom][pos] = 0
                    #maybe 1000?

    for chrom in denominator_dict:
        for pos in denominator_dict[chrom]:
            if not pos in numerator_dict[chrom]:
                ratio_dict[chrom][pos] = 0
                    #maybe 0.001?

    return ratio_dict

def write_dict_to_wig(ratio_dict, outputwig):

    outputwig_file = open(outputwig, 'w')

    for chrom in ratio_dict:
        headerline = 'variableStep chrom=%s\n'%(chrom)
        outputwig_file.write(headerline)

        for pos in sorted(ratio_dict[chrom].keys()):
            outputline = '%d\t%f\n'%(pos, ratio_dict[chrom][pos])
            outputwig_file.write(outputline)

    outputwig_file.close()


def make_avg(numerator_dict, denominator_dict):

    avg_dict = {}
    for chrom in numerator_dict:

        avg_dict[chrom] = {}
        for pos in numerator_dict[chrom]:
            if pos in denominator_dict[chrom]:
                avg_dict[chrom][pos] = (numerator_dict[chrom][pos] + denominator_dict[chrom][pos])/2
            else:
                avg_dict[chrom][pos] = 0


    for chrom in denominator_dict:
        for pos in denominator_dict[chrom]:
            if not pos in numerator_dict[chrom]:
                avg_dict[chrom][pos] = 0

    return avg_dict


def write_dict_to_maplot(ratio_dict, avg_dict, outputma):

    outputma_file = open(outputma, 'w')

    for chrom in ratio_dict:

        for pos in sorted(ratio_dict[chrom].keys()):
            outputline = '%s\t%d\t%f\t%f\n'%(chrom, pos, ratio_dict[chrom][pos], avg_dict[chrom][pos])
            outputma_file.write(outputline)


    outputma_file.close()

def create_arrays(ratio_dict, avg_dict):
    #in MA plot x values are avg reads, y values are fold-change
    x_array = []
    y_array = []
    label_array = []

    for chrom in ratio_dict:

        for pos in sorted(ratio_dict[chrom].keys()):
            x_array.append(avg_dict[chrom][pos])
            y_array.append(ratio_dict[chrom][pos])
            label_array.append(str(chrom) + '_' + str(pos))

    return x_array, y_array, label_array

def log_hover_scatter(x_array, y_array, label_array, out_prefix, xlabel='', ylabel=''):
    from bokeh.plotting import figure, output_file, save, ColumnDataSource, gridplot
    from bokeh.io import push_notebook, output_notebook, show
    from bokeh.models import Range1d
    from bokeh.models import HoverTool
    from collections import OrderedDict
    output_notebook()

    bokeh_black = (0, 0, 0)
    bokeh_orange = (230, 159, 0)
    bokeh_skyBlue = (86, 180, 233)
    bokeh_bluishGreen = (0, 158, 115)
    bokeh_yellow = (240, 228, 66)
    bokeh_blue = (0, 114, 178)
    bokeh_vermillion = (213, 94, 0)
    bokeh_reddishPurple = (204, 121, 167)

    output_file("%s.html" % (out_prefix))

    TOOLS = "pan,wheel_zoom,reset,save"

    source = ColumnDataSource(data=dict(x=x_array, y=y_array, label=label_array))

    p1 = figure(x_axis_label=xlabel, y_axis_label=ylabel,
                y_axis_type="log", x_axis_type="log", tools=TOOLS, toolbar_location="right")
    r1 = p1.circle('x', 'y', size=5, source=source, color=bokeh_black, name='main')
    p1.x_range = Range1d(start=0.00001, end=100)
    p1.y_range = Range1d(start=0.00001, end=100)
    hover = HoverTool(names=['main'])
    hover.tooltips = [('', '@label')]
    p1.add_tools(hover)
    save(p1)
    show(p1)


numerator, denominator, outputwig, outputma, out_prefix = sys.argv[1:]
    #take in arguments from command line
numerator_dict = parse_wig(numerator)
denominator_dict = parse_wig(denominator)
    #create dictionaries
ratio_dict = make_ratio(numerator_dict, denominator_dict)
    #create dicitonary for ratios
write_dict_to_wig(ratio_dict, outputwig)
    #write dictionary to a wig file output specified
avg_dict = make_avg(numerator_dict, denominator_dict)
    #create dictionary for average reads at position across the two files
write_dict_to_maplot(ratio_dict, avg_dict, outputma)
    #write values to a tab-delimeted file for MA plots
x_array, y_array, label_array = create_arrays(ratio_dict, avg_dict)
    #create arrays for bokeh plotting
log_hover_scatter (x_array, y_array, label_array, out_prefix, xlabel='', ylabel='')
    #plot MA plots in bokeh