#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AUTHOR

    Pedro Cerqueira
    github: @pedrorvc
    
DESCRIPTION


"""

import csv
import os
import argparse

import xml.etree.ElementTree as ET
from xml.dom import minidom
from matplotlib import cm
from collections import OrderedDict

from Bio import SeqIO

def listdir_fullpath(path):
    """ Gets the full path of the files from a directory

        Args:
            path (str): full path to a directory

        Returns:
            list containing the full path of every file contained in the input directory
    
    """
    
    return [os.path.join(path, f) for f in os.listdir(path)]


def prettify(elem):
    """Return a pretty-printed XML string for the Element.
    """
    rough_string = ET.tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    
    return reparsed.toprettyxml(indent="  ")


def ring_attributes(colour, name, position):
    """
    """
    
    ring_attrs = {"colour" : colour,
                 "name": name,
                 "position" : position,
                 "upperInt" : "90",
                 "lowerInt" : "70",
                 "legend" : "yes",
                 "size" : "30",
                 "labels" : "no",
                 "blastType" : "blastn"}
    
    return ring_attrs


def annotation_ring_attributes(position):
    """
    """
    
    annotation_ring_attrs = {"colour" : '172,14,225',
                 "name": 'null',
                 "position" : position,
                 "upperInt" : "70",
                 "lowerInt" : "50",
                 "legend" : "yes",
                 "size" : "30",
                 "labels" : "no",
                 "blastType" : "blastn"}
    
    return annotation_ring_attrs


def create_feature_attrs(label, colour, decoration, start, stop):
    """
    """
    
    feature_element_attrs = {'label' : label,
                             'colour' : colour,
                             'decoration' : decoration}
    
    feature_range_element_attrs = {'start' : start,
                                   'stop' : stop}
    
    return feature_element_attrs, feature_range_element_attrs


def create_annotation_ring_tsv(annotation_ring, annotation_file):
    """
    """
    
    with open(annotation_file) as tsvfile:
      reader = csv.DictReader(tsvfile, dialect='excel-tab')
      
      # Obtain the annotations from the file contents
      for row in reader:
          start = row['#START']
          stop = row['STOP']
          label = row['Label']
          colour = row['Colour']
          decoration = row['Decoration']
            
          # Create xml attributes
          feature_element_attrs, feature_range_element_attrs = create_feature_attrs(label, colour, decoration, start, stop)
          
          # Create xml elements
          feature_element = ET.SubElement(annotation_ring, 'feature', attrib=feature_element_attrs)        
          feature_range_element = ET.SubElement(feature_element, 'featureRange', attrib=feature_range_element_attrs)
          

#def create_gbk_genes_of_interest_feature_elements_concat(annotation_ring, records, genes_of_interest):
#    """
#    """
#    
#    with open(genes_of_interest, "r") as f:
#        genes = f.readlines()
#        genes = [gene.rstrip() for gene in genes]
#    
#    for seq_record in records:
##        ass_num = seq_record.id
#
#        for f in seq_record.features:
#            
#            if f.type == 'CDS':
#                            
#                if 'gene' in f.qualifiers and f.qualifiers['gene'][0] in genes:
#                    label = f.qualifiers['gene'][0]
#                elif 'product' in f.qualifiers and f.qualifiers['product'][0] in genes:
#                    product = f.qualifiers['product'][0]
#                    label = product
#                else:
#                    #print(ass_num, "No gene or product tag!\n")
#                    continue
#                    
#                
#                start = str(f.location.start.position)
#                end = str(f.location.end.position)
#                strand = f.location.strand
#                    
#                    
#                if strand == -1:
#                    decoration = 'counterclockwise-arrow'
#                
#                elif strand == 1:
#                    decoration = 'clockwise-arrow'
#                
#                # Create xml attributes
#                feature_element_attrs, feature_range_element_attrs = create_feature_attrs(label, "black", decoration, start, end)
#                  
#                # Create xml elements
#                feature_element = ET.SubElement(annotation_ring, 'feature', attrib=feature_element_attrs)        
#                feature_range_element = ET.SubElement(feature_element, 'featureRange', attrib=feature_range_element_attrs)
      
              
#def create_gbk_feature_elements_concat(annotation_ring, records):
#    """
#    """
#    
#    for seq_record in records:
##        ass_num = seq_record.id
#        
#        for fea in seq_record.features:
#            if fea.type == 'CDS':
#                start = str(fea.location.start.position)
#                end = str(fea.location.end.position)
#                strand = fea.location.strand
#                
#                if 'gene' in fea.qualifiers:
#                    label = str(fea.qualifiers['gene'][0])
#                elif 'product' in fea.qualifiers:
#                    product = fea.qualifiers['product'][0]
#                    label = str(product)
#                else:
##                    print(ass_num, "No gene or product tag!\n")
#                    continue
#                
#                if strand == -1:
#                    decoration = 'counterclockwise-arrow'
#                
#                elif strand == 1:
#                    decoration = 'clockwise-arrow'
#                    
#                # Create xml attributes
#                feature_element_attrs, feature_range_element_attrs = create_feature_attrs(label, "black", decoration, start, end)
#                  
#                # Create xml elements
#                feature_element = ET.SubElement(annotation_ring, 'feature', attrib=feature_element_attrs)        
#                feature_range_element = ET.SubElement(feature_element, 'featureRange', attrib=feature_range_element_attrs)
                
                
def annotation_ring_feature_elements_gbk_concat(annotation_ring, record, genome_size):
    """
    """
    
    if type(genome_size) == int:

    
        for fea in record.features:
            
            if fea.type == 'source':
                size = fea.location.end.position
    
          
            if fea.type == 'CDS':
                start = str(fea.location.start.position)
                end = str(fea.location.end.position)
                strand = fea.location.strand
                
                if 'gene' in fea.qualifiers:
                    label = str(fea.qualifiers['gene'][0])
                elif 'product' in fea.qualifiers:
                    product = fea.qualifiers['product'][0]
                    label = str(product)
                else:
        #                    print(ass_num, "No gene or product tag!\n")
                    continue
                
                if strand == -1:
                    decoration = 'counterclockwise-arrow'
                
                elif strand == 1:
                    decoration = 'clockwise-arrow'
                    
                # Create xml attributes
                feature_element_attrs, feature_range_element_attrs = create_feature_attrs(label, "black", decoration, start, end)
                  
                # Create xml elements
                feature_element = ET.SubElement(annotation_ring, 'feature', attrib=feature_element_attrs)        
                feature_range_element = ET.SubElement(feature_element, 'featureRange', attrib=feature_range_element_attrs)
        
        genome_size += size
            

def annotation_ring_feature_elements_genes_of_interest_gbk_concat(annotation_ring, record, genes, genome_size=False):
    """
    """
#    print(record.id)
    
    if type(genome_size) == int:
        
        for f in record.features:
            if f.type == "source":
                size = f.location.end.position
    #            print(size)
        
            if f.type == 'CDS':
                            
                if 'gene' in f.qualifiers and f.qualifiers['gene'][0] in genes:
                    label = f.qualifiers['gene'][0]
                elif 'product' in f.qualifiers and f.qualifiers['product'][0] in genes:
                    product = f.qualifiers['product'][0]
                    label = product
                else:
                    #print(ass_num, "No gene or product tag!\n")
                    continue
                    
                
                start = str(f.location.start.position + genome_size)
                end = str(f.location.end.position + genome_size)
                strand = f.location.strand
                    
                    
                if strand == -1:
                    decoration = 'counterclockwise-arrow'
                
                elif strand == 1:
                    decoration = 'clockwise-arrow'
                
#                print(label, start, end)
                
                # Create xml attributes
                feature_element_attrs, feature_range_element_attrs = create_feature_attrs(label, "black", decoration, start, end)
                  
                # Create xml elements
                feature_element = ET.SubElement(annotation_ring, 'feature', attrib=feature_element_attrs)        
                feature_range_element = ET.SubElement(feature_element, 'featureRange', attrib=feature_range_element_attrs)
                
        genome_size += size
        
    return genome_size
      

def create_annotation_ring_gbk_concat(annotation_ring, annotation_file, genes_of_interest, records):
    """
    """
    
    if genes_of_interest != []:
        
        with open(genes_of_interest, "r") as f:
            genes = f.readlines()
            genes = [gene.rstrip() for gene in genes]
        
        for seq_record in records:
            annotation_ring_feature_elements_genes_of_interest_gbk_concat(annotation_ring, seq_record, genes)

    else:
        
        for seq_record in records:
            annotation_ring_feature_elements_gbk_concat(annotation_ring, seq_record)

#def create_gbk_feature_elements_contigs(annotation_ring, records, seq_records_dict, content_dict):
#    """
#    """
#    
#    genome_size = 0
#    for i in range(1, len(records)+1):                   
#        ord_record = seq_records_dict[content_dict[str(i)]]
#        
#        for fea in ord_record.features:
#            if fea.type == 'source':
#                size = fea.location.end.position
#                
#            
#            if fea.type == 'CDS':
#                
#                start = str(fea.location.start.position + genome_size)
#                end = str(fea.location.end.position + genome_size)
#                strand = fea.location.strand
#            
#                if 'gene' in fea.qualifiers:
#                    label = fea.qualifiers['gene'][0]
#                elif 'product' in fea.qualifiers:
#                    product = fea.qualifiers['product'][0]
#                    label = product
#                else:
#    #                print(ass_num, "No gene or product tag!\n")
#                    continue
#                
#                if strand == -1:
#                    decoration = 'counterclockwise-arrow'
#                
#                elif strand == 1:
#                    decoration = 'clockwise-arrow'
#                    
#                # Create xml attributes
#                feature_element_attrs, feature_range_element_attrs = create_feature_attrs(str(label), 'black', decoration, str(start), str(end))
#                  
#                # Create xml elements
#                feature_element = ET.SubElement(annotation_ring, 'feature', attrib=feature_element_attrs)        
#                feature_range_element = ET.SubElement(feature_element, 'featureRange', attrib=feature_range_element_attrs)
#            
#        genome_size += size

#def create_gbk_feature_genes_of_interest_elements_contigs(annotation_ring, records, genes_of_interest, seq_records_dict, content_dict):
#    """
#    """
#    
#    with open(genes_of_interest, "r") as f:
#        genes = f.readlines()
#        genes = [gene.rstrip() for gene in genes]
#    
#    genome_size = 0
#    for i in range(1, len(records)+1):                   
#        ord_record = seq_records_dict[content_dict[str(i)]]
#        
#        for f in ord_record.features:
#        
#            if f.type == 'source':
#                size = f.location.end.position
#            
#            if f.type == 'CDS':
#                        
#                if 'gene' in f.qualifiers and f.qualifiers['gene'][0] in genes:
#                    label = f.qualifiers['gene'][0]
#                elif 'product' in f.qualifiers and f.qualifiers['product'][0] in genes:
#                    product = f.qualifiers['product'][0]
#                    label = product
#                else:
#                    #print(ass_num, "No gene or product tag!\n")
#                    continue
#                    
#                start = str(f.location.start.position + genome_size)
#                end = str(f.location.end.position + genome_size)
#                strand = f.location.strand
#                    
#                if strand == -1:
#                    decoration = 'counterclockwise-arrow'
#                
#                elif strand == 1:
#                    decoration = 'clockwise-arrow'
#                                                    
#                # Create xml attributes
#                feature_element_attrs, feature_range_element_attrs = create_feature_attrs(label, "black", decoration, start, end)
#                  
#                # Create xml elements
#                feature_element = ET.SubElement(annotation_ring, 'feature', attrib=feature_element_attrs)        
#                feature_range_element = ET.SubElement(feature_element, 'featureRange', attrib=feature_range_element_attrs)
#                
#        genome_size += size
        

def create_annotation_ring_gbk_contigs(annotation_ring, annotation_file, records, genes_of_interest, contig_order):
    """
    """
    
    if contig_order != []:
          
        with open(contig_order) as tsvfile:
            reader = csv.DictReader(tsvfile, dialect='excel-tab')
        
            content_dict = OrderedDict()
            for r in reader:
                content_dict[r["order"]] = r["contig"]
    

        seq_records_dict = OrderedDict()
        for record in records:            
            seq_records_dict[record.id] = record
            
        if genes_of_interest != []:
            
            with open(genes_of_interest, "r") as f:
                genes = f.readlines()
                genes = [gene.rstrip() for gene in genes]
                
            genome_size = 0
            for i in range(1, len(records)+1):                                   
                ord_record = seq_records_dict[content_dict[str(i)]]
                                
                gsize = annotation_ring_feature_elements_genes_of_interest_gbk_concat(annotation_ring, ord_record, genes, genome_size)
                
                genome_size = gsize
        else:
            
            genome_size = 0
            for i in range(1, len(records)+1):                                   
                ord_record = seq_records_dict[content_dict[str(i)]]
                
                gsize = annotation_ring_feature_elements_gbk_concat(annotation_ring, ord_record, genome_size)
                
                genome_size = gsize

    else:
    
        if genes_of_interest != []:
            
            with open(genes_of_interest, "r") as f:
                genes = f.readlines()
                genes = [gene.rstrip() for gene in genes]

            
            for seq_record in records:
                annotation_ring_feature_elements_genes_of_interest_gbk_concat(annotation_ring, seq_record, genes)
    
        else:
            for seq_record in records:
                annotation_ring_feature_elements_gbk_concat(annotation_ring, seq_record)


def write_xml(root_elem, output_file):
    """
    """
    #brig_params.xml

    xml_file = ET.tostring(root_elem, encoding='utf8').decode('utf8')
    
    pretty_xml_file = minidom.parseString(xml_file).toprettyxml(indent='    ')
    
    with open(output_file, "w") as f:
        f.write(pretty_xml_file)




# Get information 

## Create root
### Root attributes (BRIG)
        
def create_root_element(blast_options, legend_position, query_file, 
                               output_folder, image_output_file, title, image_format):
    """
    """

    root_attrs = {"blastOptions" : blast_options,
                  "legendPosition" : legend_position,
                  "queryFile" : query_file,
                  "outputFolder" : output_folder,
                  "blastPlus" : "yes",
                  "outputFile" : os.path.join(output_folder, image_output_file),
                  "title" : title,
                  "imageFormat" : image_format,
                  "queryFastaFile" : query_file,
                  "cgXML" : os.path.join(output_folder + "/scratch", os.path.basename(query_file) + ".xml")}
    
    
    root = ET.Element('BRIG', attrib=root_attrs)
    
    return root

#### Create root children

# Create cgview_settings
    
def create_cgview_settings_element(root, height, width):
    """
    """

    cgview_settings_attrs = {"arrowheadLength" : "medium",
                             "backboneColor" : "black",
                             "backboneRadius" : "600",
                             "backboneThickness" : "medium",
                             "backgroundColor" : "white",
                             "borderColor" : "black",
                             "featureSlotSpacing" : "medium",
                             "featureThickness" : "30",
                             "giveFeaturePositions" : "false",
                             "globalLabel" : "true",
                             "height" : height,
                             "isLinear" : "false",
                             "labelFont" : "SansSerif,plain,25",
                             "labelLineLength" : "medium",
                             "labelLineThickness" : "medium",
                             "labelPlacementQuality" : "best",
                             "labelsToKeep" : "1000",
                             "longTickColor" : "black",
                             "minimumFeatureLength" : "medium",
                             "moveInnerLabelsToOuter" :"true",
                             "origin" : "12",
                             "rulerFont" : "SansSerif,plain,35",
                             "rulerFontColor" : "black",
                             "rulerPadding" : "40",
                             "rulerUnits" : "bases",
                             "shortTickColor" : "black",
                             "shortTickThickness" : "medium",
                             "showBorder" : "false",
                             "showShading" : "true",
                             "showWarning" : "false",
                             "tickDensity" : "0.2333",
                             "tickThickness" : "medium",
                             "titleFont" : "SansSerif,plain,45",
                             "titleFontColor" : "black",
                             "useColoredLabelBackgrounds" : "false",
                             "useInnerLabels" : "true",
                             "warningFont" : "Default,plain,35",
                             "warningFontColor" : "black",
                             "width" : width,
                             "zeroTickColor" : "black",
                             "tickLength" : "medium"}


    cgview_settings = ET.SubElement(root, 'cgview_settings', attrib=cgview_settings_attrs)
    
    
#    return cgview_settings

# Create brig_settings
    
def create_brig_settings_element(root, java_memory):
    """
    """

    brig_settings_attrs = {"Ring1" : "172,14,225",
                           "Ring2" : "222,149,220",
                           "Ring3" : "161,221,231",
                           "Ring4" : "49,34,221",
                           "Ring5" : "116,152,226",
                           "Ring6" : "224,206,38",
                           "Ring7" : "40,191,140",
                           "Ring8" : "158,223,139",
                           "Ring9" : "226,38,122",
                           "Ring10" :"211,41,77",
                           "defaultUpper" : "70",
                           "defaultLower" : "50",
                           "defaultMinimum" : "50",
                           "genbankFiles" : "gbk,gb,genbank",
                           "fastaFiles" : "fna,faa,fas,fasta,fa",
                           "emblFiles" : "embl",
                           "blastLocation" : "",
                           "divider" : "3",
                           "multiplier" : "3",
                           "memory" : java_memory,
                           "defaultSpacer" : "0"}
    
    brig_settings = ET.SubElement(root, 
                                  "brig_settings", 
                                  attrib=brig_settings_attrs)


## Create special

def create_special_element(root):
    """
    """
    gc_content_special = ET.SubElement(root, 'special', attrib={'value' : 'GC Content'})
    gc_skew_special = ET.SubElement(root, 'special', attrib={'value' : 'GC Skew'})

# Create reference dir element

def create_reference_directory_element(root, reference_directory):
    """
    """

    ref_dir = ET.SubElement(root, 
                            "refDir", 
                            attrib={"location" : reference_directory})
    
    ref_dir_list = listdir_fullpath(reference_directory)
    
    for f in ref_dir_list:
        ref_file = ET.SubElement(ref_dir, 
                                 "refFile", 
                                 attrib={"location" : f})





def create_annotation_ring(root, reference_directory, annotation_file, genes_of_interest, contig_order):
    """
    """
    
    ring_position = len(os.listdir(reference_directory)) + 2
    
    annotation_ring = ET.SubElement(root, 'ring', attrib=annotation_ring_attributes(str(ring_position)))
    
    if list(SeqIO.parse(annotation_file, "genbank")) == []:
        create_annotation_ring_tsv(annotation_ring, annotation_file)
    
    else:
        
        records = [r for r in SeqIO.parse(annotation_file, "genbank")]
        
        if len(records) > 1:
            create_annotation_ring_gbk_contigs(annotation_ring, annotation_file, records, genes_of_interest, contig_order)
        else:
            create_annotation_ring_gbk_concat(annotation_ring, annotation_file, genes_of_interest, records)
            


## Create remaining rings
        
def create_ring_element(root, reference_directory, colormap):
    """
    """

#    reference_directory = "/home/pcerqueira/Lab_Software/misc_scripts/reference_dir"
    #    reference_directory = "/home/pcerqueira/Lab_Software/misc_scripts/ordered_reference_dir"

    
    ref_dir_list = listdir_fullpath(reference_directory)    
    
    # Gets the colormap from matplotlib with as many colors as the number of files
    cmap = cm.get_cmap(colormap, len(ref_dir_list))
    # cmap = cm.get_cmap('viridis', len(ref_dir_list))
    list_colormap = cmap.colors.tolist()
    
    # Remove the fourth element (transparency) because it is not necessary
    colors_to_use = []
    for l in list_colormap:
        convert = [round(x * 255) for x in l]
        convert.pop()
        colors_to_use.append(convert)
    
    #reversed_colors_to_use = colors_to_use[::-1]
    
    # Check if the ueser provided ao order for the rings
    has_digit = [os.path.basename(x).split("_")[0].isdigit() for x in ref_dir_list] 
    
    if True in has_digit:
        ring_positions = [os.path.basename(x).split("_")[0] for x in ref_dir_list]
        ring_positions.sort(reverse=True)
        ref_dir_list.sort(reverse=True)
        
        for ring in range(len(ref_dir_list)):
            
            ring_position = int(ring_positions[ring]) + 1
            
            ring_color = ",".join([str(e) for e in colors_to_use[ring]])
            
            ring_name = os.path.basename(ref_dir_list[ring]).split("_")[1]
                        
            ring_number_element = ET.SubElement(root, 
                                    'ring',
                                    ring_attributes(ring_color, ring_name, str(ring_position)))
            
            ring_sequence_element = ET.SubElement(ring_number_element, 
                                      "sequence", 
                                      attrib={"location" : ref_dir_list[ring]})
        
        
    else:
        # Sort files by lowercase
        ref_dir_list.sort(key=lambda y: y.lower())
 
        # The number of rings starts at 2 due to the GC Content and GC Skew
        ring_number = len(ref_dir_list) + 1
        for ring in range(len(ref_dir_list)):
            
            ring_color = ",".join([str(e) for e in colors_to_use[ring]])
            
            ring_name = os.path.basename(ref_dir_list[ring]).split("_")[0]
            
            ring_number_element = ET.SubElement(root, 
                                                'ring',
                                                ring_attributes(ring_color, ring_name, str(ring_number)))
            
            ring_sequence_element = ET.SubElement(ring_number_element, 
                                                  "sequence", 
                                                  attrib={"location" : ref_dir_list[ring]})
            
            ring_number -= 1
    
  
## Create special rings

def create_special_ring_element(root):
    """
    """
    
    # Create ring attributes
    gc_content_ring_attrs = ring_attributes('225,0,0', "GC Content", "0")
    gc_skew_ring_attrs = ring_attributes('225,0,0', "GC Skew", "1")
    
    # Add ring element to root
    gc_skew_ring = ET.SubElement(root, 'ring', attrib=gc_skew_ring_attrs)
    gc_content_ring = ET.SubElement(root, 'ring', attrib=gc_content_ring_attrs)
    
    # Add sequence element to ring
    gc_content_location = ET.SubElement(gc_content_ring, 'sequence', attrib={'location' : 'GC Content'})
    gc_skew_location = ET.SubElement(gc_skew_ring, 'sequence', attrib={'location' : 'GC Skew'})





def main(query_file, reference_directory, output_folder, output_file, image_output_file, title, annotation_file,
         genes_of_interest, contig_order, blast_options, legend_position, image_format, height, width, java_memory, colormap):
    """
    """
    
    # Create root element
    root = create_root_element(blast_options, legend_position, query_file, 
                               output_folder, image_output_file, title, image_format)
    
#    root = create_root_element(query_file, output_folder, image_output_file, title,
#                               blast_options, legend_position, image_format)
    
    create_cgview_settings_element(root, height, width)
    
    create_brig_settings_element(root, java_memory)
    
    create_special_element(root)

    create_reference_directory_element(root, reference_directory)
        
    if annotation_file:
        create_annotation_ring(root, reference_directory, annotation_file, genes_of_interest, contig_order)
    
    create_ring_element(root, reference_directory, colormap)
        
    create_special_ring_element(root)
    
    write_xml(root, output_file)

def parse_arguments():
    
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('-q', '--query', type=str, required=True, dest='query_file',
                        help='Path to the query/reference FASTA file.')
    
    parser.add_argument('-rfd', '--ref_dir', type=str, required=True, dest='reference_directory',
                        help='Path to the directory where the FASTA file to compare against the reference.')
    
    parser.add_argument('-od', '--out_dir', type=str, required=True, dest='output_folder',
                        help='Path to the output directory for the results of BRIG.')
    
    parser.add_argument('-of', '--out_file', type=str, required=True, dest='output_file',
                        help='Path to the output of this script.')

    parser.add_argument('-oi', '--out_img', type=str, required=True, dest='image_output_file',
                        help='Path to the output file of the resulting image of BRIG.')

    parser.add_argument('-t', '--title', type=str, required=True, dest='title',
                        help='Title of the resulting image from BRIG.')
    
    parser.add_argument('-a', '--annotation', type=str, required=False, dest='annotation_file', default=False,
                        help='File containing annotations for the reference genome. '
                        'The annoation file can be a tab-delimited file (.tsv) or a Genbank format file (.gbk, .gb)')
    
    parser.add_argument('--genes', type=str, required=False, dest='genes_of_interest', default=[],
                        help='File containing a list of specific genes (one gene per line) to search when a Genbank annotation file is provided. ')
   
    parser.add_argument('--contig-order', type=str, required=False, dest='contig_order', default=[],
                        help='Tab-delimited file containing the order of the contigs when a Genbank (divided by contigs) annotation file is provided. '
                             'Example:  order    contig '
                                          '1     Contig8')

    parser.add_argument('-b', '--blast_options', type=str, required=False, dest="blast_options", default="-evalue 0.001 -num_threads 6",
                        help='Options for running BLAST.')
    
    parser.add_argument('-l', '--legend_pos', type=str, required=False, dest="legend_position", default="middle-right",
                        help='Positon of the legend on the resulting image.'
                        'The options available are upper, center or lower, '
                        'paired with left, center or right')
    
    parser.add_argument('-if', '--image_format', type=str, required=False, dest="image_format", default="jpg",
                        help='Format of the resulting image file.'
                        'The available options are: jpg, png, svg or svgz.')

    parser.add_argument('-ht', '--height', type=str, required=False, dest="height", default="3000",
                        help='Height (in pixels) of the resulting image.')
    
    parser.add_argument('-wd', '--width', type=str, required=False, dest="width", default="3000",
                        help='Width (in pixels) of the resulting image.')

    parser.add_argument('-jm', '--java_memory', type=str, required=False, dest="java_memory", default="1500",
                        help='Amount of memory (in bytes) that Java is allowed to use for BRIG.')

    parser.add_argument('-cm', '--colormap', type=str, required=False, dest="colormap", default="viridis",
                        help='Colormap from matplotlib to use for the color of the rings. '
                        'The available options are: viridis, plasma, inferno, magma and cividis.'
                        'More options for colormaps at: '
                        'https://matplotlib.org/users/colormaps.html')

    args = parser.parse_args()
    
    return [args.query_file, args.reference_directory, args.output_folder, args.output_file,
            args.image_output_file, args.title, args.annotation_file, args.genes_of_interest, args.contig_order,
            args.blast_options, args.legend_position, args.image_format, args.height, args.width, args.java_memory, args.colormap]
    
if __name__ == '__main__':
    
    args = parse_arguments()
    main(args[0], args[1], args[2], args[3], args[4], args[5], args[6],
         args[7], args[8], args[9], args[10], args[11], args[12], args[13],
         args[14], args[15])
    
