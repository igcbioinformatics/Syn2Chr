#!/nfs/production/mousegenomes/bin/python-3.5.6/python
# Version 2018-4-22: 
#    Repeat filtering added. In the same 10K adjacent windows, if there are more than 5 BLAST hits going to more than 3 different chromosome, all hits will be removed.
# Version 2018-4-24:
#    Change output format of hit_list
#    Argument processing changed to argparse library
#    offset of scaffold added, offer better alignment between chromosomes
# Version 2018-10-5:
#   Restructure according to new pipeline. Now Scaffold files BLAST onto Repeatmarsked GRCm. Rename from crosslink.py to synteny.py
# Version 2018-10-15:
#   Update SvgObj.ChrBox module, add option for other species, and complex loci highlighting
# Version 2018-12-01:
#   Add chromosome G-banding prediction for scaffold
# Todo:
#    generate cmap and xmap for irysview
#    internal .fasta.fai
#    Should support referene genome other than Mus Mus

import os
import sys
import argparse
import img2pdf
import operator
import statistics
from wand.image import Image

########################################## Functions #############################################
class SvgObj:
    def __init__(self): # colours available for up to 40 different chromosome's.
        self.chr_color = {'chr1':'#ffaaaa','chr2':'#ff9955','chr3':'#ffd42a','chr4':'#ccff00','chr5':'#66ff00','chr6':'#00ffcc','chr7':'#000080','chr8':'#80b3ff','chr9':'#b380ff','chr10':'#f4d7ee','chr11':'#ff0000','chr12':'#803300','chr13':'#a0892c','chr14':'#677821','chr15':'#5aa02c','chr16':'#165016','chr17':'#0066ff','chr18':'#214478','chr19':'#c837ab','chr20':'#ef9504', 'chr21':'#4c2418','chr22':'#e8e8a0','chr23':'#5e5e31','chr24':'#315e53','chr25':'#79fce6','chr26':'#f442e5','chr27':'#751d6d','chr28':'#a3609d','chr29':'#c3bee2','chr30':'#919607','chr31':'#d5eaa4','chr32':'#d69f77','chr33':'#28ff89','chr34':'#fdff87','chr35':'#a00303','chr36':'#bfaddb','chr37':'#aa874e','chr38':'#e28cba','chrx':'#4400aa','other':'#b0b0b0'}

    def scale(self,x,y,xscl,xlen, font_size = 25):
        # Generate Scale line in SVG format. 
        # x,y is starting point of scale line. xscl is length of the line, xlen is the length of genome/scaffold
        zoom, m_scale, Svg_line, Svg_text, text = 0, 0, '', '', ''
        # adjust ruler scale and round up, variable 'zoom' is bp(xlen) per ruler line
        zoom = xlen/10 if xscl > 1000 else xlen*100/xscl
        m_scale = {'1':10,'2':2,'3':2,'4':5,'5':5,'6':5,'7':5,'8':10,'9':10}[str(zoom)[0]]
        zoom = 10**int(len(str(zoom))-1) * m_scale 
        for cx in range(xlen/zoom):
            Svg_line = Svg_line + 'M' + str(x + xscl*cx*zoom/xlen) + ' ' + str(y) + 'v 100 '
            for dx in range(1,m_scale): Svg_line = Svg_line + 'M' + str(x + xscl*cx*zoom/xlen + xscl*dx*zoom/xlen/m_scale) + ' ' + str(y) + 'v 50 '
            
            if len(str(zoom)) >= 7: text = str(int(zoom * cx /1e6)) + 'M' 
            elif len(str(zoom)) >= 4: text = str(int(zoom * cx /1e3)) + 'K'
            else: text = str(zoom * cx)

            Svg_text = Svg_text + '<text x="' + str(x + xscl*cx*zoom/xlen-20) +'" y="' + str(y+190) +'" fill="black" font-size="'+ str(font_size) + '">' + text +'</text>\n'
        return '<path d= "M 1 400' + Svg_line +' " fill = "none" stroke = "black" stroke-width="1" />\n' + Svg_text
    
    @classmethod
    def ColorGradient(self,count,maxres,ctype = 'HM'):
        Scheme = {
            'BW': [[0,255,255],0,[0,0,0],50,[128,128,128],100,[255,255,255],101,[255,0,255],[255,0,255]], #black-white
            'C1': [[255,255,255],3,[255,255,255],50,[255,0,0],100,[0,0,0],101,[0,0,0],[0,0,0]], #white-red-black
            'C2': [[255,255,200],10,[255,255,200],60,[0,255,0],100,[0,0,0],101,[0,0,0],[0,0,0]], #yellow-blue-black
            'HM': [[0,0,0],0,[0,0,0],20,[200,0,0],70,[200,200,100],100,[240,240,255],[255,0,255]] # heat map, black-red-yellow-white
            } [ctype]
        percentage = int(count*100/maxres)
        if percentage < Scheme[1]: return('rgb('+','.join(str(cx) for cx in Scheme[0])+')')
        if percentage > Scheme[7]: return('rgb('+','.join(str(cx) for cx in Scheme[9])+')')
        elif percentage <= Scheme[3]:
            return('rgb('+str(int(Scheme[2][0]+percentage*(Scheme[4][0]-Scheme[2][0])/(Scheme[3]-Scheme[1])))+','+str(int(Scheme[2][1]+percentage*(Scheme[4][1]-Scheme[2][1])/(Scheme[3]-Scheme[1])))+','+str(int(Scheme[2][2]+percentage*(Scheme[4][2]-Scheme[2][2])/(Scheme[3]-Scheme[1])))+')')
        elif percentage >= Scheme[5]:
            return('rgb('+str(int(Scheme[6][0]+(percentage-Scheme[5])*(Scheme[8][0]-Scheme[6][0])/(Scheme[7]-Scheme[5])))+','+str(int(Scheme[6][1]+(percentage-Scheme[5])*(Scheme[8][1]-Scheme[6][1])/(Scheme[7]-Scheme[5])))+','+str(int(Scheme[6][2]+(percentage-Scheme[5])*(Scheme[8][2]-Scheme[6][2])/(Scheme[7]-Scheme[5])))+')')
        else:
            return('rgb('+str(int(Scheme[4][0]+(percentage-Scheme[3])*(Scheme[6][0]-Scheme[4][0])/(Scheme[5]-Scheme[3])))+','+str(int(Scheme[4][1]+(percentage-Scheme[3])*(Scheme[6][1]-Scheme[4][1])/(Scheme[5]-Scheme[3])))+','+str(int(Scheme[4][2]+(percentage-Scheme[3])*(Scheme[6][2]-Scheme[4][2])/(Scheme[5]-Scheme[3])))+')')
    def head(self, width, height): #both value in pixel
        return '<?xml version="1.0" standalone="no"?>\n<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN"\n"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">\n<svg width="' + str(width) + 'px" height="' + str(height) + '" viewBox="-120 0 ' + str(int(width+125)) + ' ' +str(height) +'" xmlns="http://www.w3.org/2000/svg" version="1.1">\n  <title>Chromosome chain</title>\n  <desc>Generaged from BLAST</desc>\n  <rect x="0" y="0" width="' + str(width) + '" height="' + str(height) + '" style="fill:white"/>'

    def ChrBox(self, x, y, text, genome, Multiref, optionP, chr_detail = []):
        # defines complex regions on ref
        complexR = ''
        chr_ = text.lower() if text.lower() in self.chr_color.keys() else 'other'
        # defines position of complex region of each chr
        if not chr_detail:
            if genome == 'MM10': # indicator of complex region on available for mus musculus
                chr_detail = {
                'MM10':{
                'chr1':[195471971,[[9400000,9600000],[21200000,21400000],[46200000,46800000],[53200000,53800000],[58600000,59000000],[82800000,83000000],[84800000,85800000],[88000000,88400000],[92400000,92600000],[93800000,94200000],[103400000,103600000],[107000000,107400000],[117600000,118200000],[139600000,140200000],[153200000,153600000],[166400000,166800000],[171400000,172000000],[173200000,173400000],[177800000,178200000],[195200000,195400000]]],
                'chr2':[182113224,[[27400000,27600000],[37200000,37400000],[77800000,78000000],[88000000,88200000],[88200000,88400000],[116400000,116600000],[122400000,122600000],[149000000,149200000],[150000000,150600000],[151000000,151600000],[154000000,154400000],[173800000,174000000],[175000000,175200000],[177000000,178200000],[181800000,182000000]]],
                'chr3':[160039680,[[3000000,3400000],[14600000,14800000],[15200000,16000000],[28800000,29000000],[31400000,31800000],[35600000,35800000],[40800000,41000000],[64000000,64800000],[77200000,77400000],[88600000,89000000],[90600000,92400000],[93600000,94800000],[97000000,97600000],[98400000,98800000],[106000000,106400000],[113200000,113400000],[142400000,142800000],[144600000,145000000]]],
                'chr4':[156508116,[[3000000,3800000],[12000000,12200000],[49200000,49600000],[60000000,60200000],[61400000,62200000],[63200000,63400000],[88600000,89000000],[112600000,114400000],[115200000,115600000],[118600000,119000000],[121200000,122800000],[137400000,137600000],[143400000,144400000],[145400000,148000000],[149600000,149800000],[156200000,156400000]]],
                'chr5':[151834684,[[3000000,4400000],[10800000,12000000],[13600000,14000000],[14800000,15800000],[23200000,23400000],[25800000,26400000],[27800000,28000000],[31200000,31600000],[73000000,73200000],[77400000,77800000],[79400000,80000000],[86800000,87400000],[88000000,88400000],[103800000,105000000],[108800000,110200000],[120800000,121000000],[137800000,138600000],[145000000,146600000],[151400000,151800000]]],
                'chr6':[149736546,[[3000000,3600000],[41800000,42200000],[43000000,43200000],[47000000,47200000],[48200000,48400000],[57000000,58800000],[59400000,60600000],[66600000,66800000],[67400000,70800000],[82600000,82800000],[85600000,86000000],[89000000,90200000],[116400000,116600000],[121400000,122400000],[123000000,124200000],[128600000,129200000],[130000000,134000000],[136200000,137000000],[137800000,138000000]]],
                'chr7':[145441459,[[3000000,3200000],[3600000,8600000],[11000000,11200000],[11600000,11800000],[12000000,12200000],[12600000,12800000],[13600000,16000000],[17200000,17400000],[18200000,18800000],[19800000,20200000],[23200000,24800000],[25400000,27200000],[27600000,28000000],[28800000,29000000],[30800000,34200000],[38200000,39600000],[41000000,44200000],[46600000,48200000],[60200000,61800000],[64000000,64200000],[84800000,86200000],[103000000,105000000],[106200000,106800000],[107800000,109000000],[131000000,131200000],[134800000,135000000],[140200000,140800000],[142200000,142400000],[143600000,143800000]]],
                'chr8':[129401213,[[3000000,3200000],[13800000,14000000],[19200000,21800000],[27200000,27400000],[28200000,28600000],[37000000,37200000],[55200000,55800000],[59800000,60000000],[69400000,69800000],[70800000,71200000],[71600000,72000000],[99600000,99800000]]],
                'chr9':[124595110,[[3000000,3400000],[35800000,36000000],[38400000,39400000],[40000000,40200000],[55800000,56200000],[78200000,78400000],[88400000,89400000],[98800000,99000000],[109400000,109800000],[114000000,114400000],[115000000,115200000],[124200000,124600000]]],
                'chr10':[130694993,[[3200000,3400000],[7200000,7600000],[22000000,22600000],[24000000,24600000],[33400000,34000000],[51400000,51600000],[55000000,55200000],[57400000,57600000],[58200000,58400000],[79000000,79600000],[82200000,82600000],[85600000,86000000],[86400000,86800000],[112000000,112200000],[118000000,118400000],[122000000,122200000]]],
                'chr11':[122082543,[[3000000,3200000],[21400000,21600000],[46600000,46800000],[48800000,49600000],[51000000,51200000],[52000000,52600000],[60400000,60800000],[71000000,71400000],[73400000,74200000],[82800000,83800000],[98600000,98800000],[99600000,99800000],[116600000,116800000],[121800000,122000000]]],
                'chr12':[120129022,[[3000000,3200000],[8400000,8600000],[17800000,21200000],[21600000,24600000],[87400000,88400000],[102600000,102800000],[103600000,104400000],[105000000,105600000],[113200000,116000000]]],
                'chr13':[120421639,[[3000000,3600000],[12600000,14000000],[19200000,19400000],[33200000,33600000],[49800000,51000000],[60800000,62800000],[65200000,68600000],[74400000,74600000],[100200000,100600000],[113000000,113200000],[119800000,120400000]]],
                'chr14':[124902244,[[4000000,4400000],[5800000,7000000],[13000000,13200000],[19000000,19800000],[26200000,26800000],[36000000,36200000],[41200000,42600000],[43000000,45000000],[51000000,52000000],[55800000,56200000],[57000000,57200000],[59200000,59400000]]],
                'chr15':[104043685,[[3000000,3200000],[8800000,9400000],[14800000,15000000],[63800000,64000000],[75000000,75600000],[77400000,77800000],[82400000,82800000],[91200000,91600000],[98200000,98400000],[101600000,101800000],[103600000,104000000]]],
                'chr16':[98207768,[[3000000,3800000],[19000000,19400000],[20800000,21200000],[32600000,32800000],[44600000,45000000],[59200000,59400000],[93400000,93800000],[98000000,98200000]]],
                'chr17':[94987271,[[3000000,3200000],[6200000,7600000],[7800000,8200000],[13000000,14200000],[17000000,17400000],[18000000,20000000],[21000000,21200000],[21600000,24000000],[27200000,27400000],[30200000,31400000],[33000000,33200000],[33800000,34600000],[35200000,37000000],[38400000,39000000],[48000000,48400000],[53600000,54000000],[71400000,71600000],[94600000,94800000]]],
                'chr18':[90702639,[[3000000,3200000],[20600000,20800000],[36800000,37000000]]],
                'chr19':[61431566,[[3000000,3400000],[8000000,8600000],[9200000,9800000],[12000000,12400000],[13000000,14000000],[33400000,33600000],[34400000,34800000],[39400000,40200000],[61000000,61400000]]],
                'chrx':[171031299,[[3200000,3400000],[4400000,4600000],[5000000,6200000],[8200000,9600000],[31000000,31400000],[32600000,32800000],[37400000,37800000],[39800000,40000000],[46000000,46200000],[50400000,50800000],[53400000,54400000],[55200000,56200000],[62400000,62600000],[70200000,70600000],[73000000,73400000],[74400000,75000000],[75600000,75800000],[77800000,78400000],[86200000,86400000],[91200000,92200000],[94600000,95000000],[102600000,103400000],[123000000,123400000],[124000000,124200000],[135000000,135600000],[136200000,137000000],[147200000,148000000],[148800000,149000000],[149600000,150000000],[153000000,153200000],[154800000,155200000],[161000000,161200000]]],
                'other':[91744698,[]]}}[genome][chr_]
                for dx in chr_detail[1]: complexR = complexR + '<rect x="' + str(x + int(dx[0]/1e5)) + '" y="' + str(y-1) + '" width="' + str(int((dx[1]-dx[0])/1e5)+1) +'" height="' + str(101) + '" style="fill:black;"/>\n'
            else:
                RefGenome = open(genome,'r')
                chrom_length_dict, gen_dict = {}, {}
                for line in RefGenome:
                    splitted_line = line.split('\t')
                    chrom_length_dict[splitted_line[0].lower()] = splitted_line[1]
                for chrom in chrom_length_dict.keys(): gen_dict[chrom] = [chrom_length_dict[chrom]]
                chr_detail = {args.Ref_genome: gen_dict}[genome][chr_]
                complexR = ''

        width = int(chr_detail[0])/1e5+1
        box_print = '<rect x="' + str(x) + '" y="' + str(y) + '" width="' + str(width) +'" height="' + str(100) + '" style="fill:'+self.chr_color[chr_]+'"/>\n'
        if text[0:3].lower() == 'chr':
            font_size = 18
            yco = int(y + 75)
            xco = int(x + width/2)
            text_print = '<text text-anchor="middle" x ="' + str(xco) + '" y="' + str(yco) + '" fill="black" font-size="' + str(font_size) +'" > ' + text + '</text>' # scaffold chr
            if x == 0:
                text_print = '<text x ="' + str(-110) + '" y="' + str(365) + '" fill="black" font-size="' + str(35) +'" > ' + text + '</text>' # ref chr
                if Multiref:
                    text_print = '<text x ="' + str(-110) + '" y="' + str(int(y + 55)) + '" fill="black" font-size="' + str(35) +'" > ' + text + '</text>' # multiple ref chr (-s option)

        else:
            font_size = 30
            yco = int(y + 45)
            xco = int(x + width/2)
            if optionP != False: 
                font_size = 18
                if optionP % 2 == 0: yco -= 15
            text_print = '<text text-anchor="middle" x ="' + str(xco) + '" y="' + str(yco) + '" fill="black" font-size="' + str(font_size) +'" > ' + text + '</text>'
        return box_print + complexR, text_print

class ChrBand:
    def __init__(self):
        self.band = { 
        # based on MM10, resolution 1M
        'chr1': '32233333300011111111|00333333333000000000|00022222222222102221|00000111111111100000|11111110000000000000|11103333333310333333|33333300122222222221|00000002333313333100|11110011110000001111|1100000011111111',
        'chr2': '32233333333333001111|11111111100000000000|00000000233333333333|23333333100022222222|20000000133333333333|31000111111110002222|22220000000011100000|03333330333310000001|11100000000000111111|11',
        'chr3': '32233333333333333310|12222222222222210000|00000003333333332000|01111111110013333333|33330000000001111100|00000011000000033333|33333330022222222210|00001111000000111111',
        'chr4': '32233333333333300133|33333333102222222222|22210000000011100000|00011111110033333333|33330000002222222100|00000122222000000011|10000000000111000000|01333333000000000',
        'chr5': '32233333333333301222|22222100000000110000|00000000000111111110|03333333333201222200|01333333333100111111|00111111000000000000|00000110111000000000|000000111111',
        'chr6': '32233333333333332000|01333331233333100000|01221022222000000000|00133333333333100222|22222200000000110000|00023333213333333000|00011000000022222221|0001111111',
        'chr7': '32233333333333300011|11111111000000111000|00000000111111100000|13333333333100000222|10000000003333333333|00111111111000000011|11100000000000000111|100000',
        'chr8': '32233333333333330111|00000000001111000000|00002222102222220001|33333330333333200000|11111000001000023333|33301111111000000000|0000111111',
        'chr9': '32233333333333100001|11110000000000000022|22220012222222200000|00011111110000000011|11101333333233333333|32000000111000000000|1111',
        'chr10':'32333333333330000011|11110000000001333333|32000000233333330000|00003333233333300000|00000000011111110003|33333333333300000000|00000111111',
        'chr11':'32233333333330000333|33000133331000003333|33310000110000000000|11100000000111000000|00333333331000000000|00122222221000000000|00',
        'chr12':'32233333333333333200|02222200000011111110|00001000000033333333|33333300000333333333|30000111000000013333|33333300000000222222',
        'chr13':'32233333333333331000|02222222210001111111|11000111111110000001|11000000001111111110|03333333333333200000|00000001111000000111|1',
        'chr14':'32233333333333300001|11111111110000000000|00013333333300022222|00000000011110000000|00000222233333333330|33333331000333333333|300',
        'chr15':'32233333333333331000|00001111110033333333|33300222220001333333|33333310022222222100|00001110000000002222|2011',
        'chr16':'32233333333333311111|10000011111100000011|11111000000000222202|22222220001333333332|333333333332000000',
        'chr17':'32233333333333001000|00222222222100000000|10000022222222220000|33333333233331000011|111000000111111',
        'chr18':'32233333333333333332|00000000012222210333|33333200013333233333|20000000111111100000|00001111111',
        'chr19':'32233333333333332000|00000122222222200022|22222221000011111110|000',
        'chrx' :'32233333333333330011|10000000222221122222|00000000222222221000|00022222210222222100|00233333333033333333|30000000003333333333|03333333333333300000|00111111100000001111|111100000000'
        }
    def BandEmulation(self,scaffold_length,window_size,group_list): # group_list: {chr:[[Ref_start,Ref_end,Que_start,Que_end,count=xxx][Ref_start,Ref_end,Que_start,Que_end,count=xxx]],chr:[...]}
        band_list, scaffold_band, scaffold_overlap, scaffold_pattern = [], ['_']*int(scaffold_length/1e6+1), [0]*int(scaffold_length/1e6+1), [0]*int(scaffold_length/1e6+1)
        for chr_ in group_list.keys(): 
            band_list += [[int(dx[2]/1e6), int(dx[3]/1e6),self.band[chr_].replace('|','')[int(min(dx[0],dx[1])/1e6):int(max(dx[0],dx[1])/1e6)]] for dx in group_list[chr_]]
        band_list = sorted(map(lambda dx: dx if dx[0]<dx[1] else [dx[1],dx[0],dx[2][::-1]],band_list),key = lambda ex: ex[0])
        for dx in band_list:
            for cx in range(dx[0],dx[1]+1): scaffold_overlap[cx] += 1
            if dx[1]-dx[0]+1 == len(dx[2]): scaffold_band[dx[0]:dx[1]+1] = list(dx[2])
            elif dx[1]-dx[0]+1 > len(dx[2]) and len(dx[2]) != 0: scaffold_band[dx[0]:dx[0]+len(dx[2])] = list(dx[2])
            elif dx[1]-dx[0]+1 < len(dx[2]): scaffold_band[dx[0]:dx[1]+1] = list(dx[2])[0:dx[1]-dx[0]+1]
        for cx in range(len(scaffold_band)):
            for dx in range(cx-window_size//2, cx+window_size//2+1):
                if dx >= 0 and dx < len(scaffold_band):
                    scaffold_pattern[cx] = scaffold_pattern[cx] + int(scaffold_band[dx]) if scaffold_band[dx] != '_' else scaffold_pattern[cx]
        return [scaffold_pattern,scaffold_overlap]

    def BandSVG(self,scaffold_length,window_size,group_list,x=0,y=0,height=100):
        output,maxres = '',window_size*3
        scaffold_pattern, scaffold_overlap = self.BandEmulation(scaffold_length,window_size,group_list)
        for cx in range(len(scaffold_pattern)):
            output += '<rect x="' + str(x+cx*10) + '" y="' + str(y) + '" width="10" height="' + str(height) + '" style="fill:'+SvgObj.ColorGradient(maxres-int(scaffold_pattern[cx]),maxres,'BW')+';"/>\n'
            output += '<rect x="' + str(x+cx*10) + '" y="' + str(y+height) + '" width="10" height="' + str(scaffold_overlap[cx]*30+10) + '" style="fill:blue"/>\n'
        return output

class IrysObj:
    def cmapHead(self, scaffold_count):
        return '# CMAP File Version:\t0.1\n# Label Channels:\t1\n# Nickase Recognition Site 1:\tGCTCTTC\n# Number of Consensus Maps:\t'+ str(scaffold_count) +'\n#h CMapId\tContigLength\tNumSites\tSiteID\tLabelChannel\tPosition\tStdDev\tCoverage\tOccurrence\n#f int\tfloat\tint\tint\tint\tfloat\tfloat\tint\tint\n'
    def xmapHead(self, R_file, Q_file):
        return '#XMAP File Version:\t0.1\n# Label Channels:\t1\n# Reference Maps From:\t./'+R_file+ '# Query Maps From:\t' + Q_file+ '\n#h XmapEntryID\tQryContigID\tRefContigID\tQryStartPos\tQryEndPos\tRefStartPos\tRefEndPos\tOrientation\tConfidence\tHitEnum\tQryLen\tRefLen\t LabelChannel\tAlignment\n#f int\tint\tint\tfloat\tfloat\tfloat\tfloat\tstring\tfloat\tstring\tfloat\tfloat\tint\tstring\n'
    def cmapRef(self,data): # data: {scaffold_name:[length,hit,hit,hit...],scaffold_name:[length,hit,hit,hit....]}
        output = ''
        for dx in data.keys():
            for cx, hit in enumerate(sorted(data[dx][1:])): output = output + '\t'.join([dx,data[dx][0],len(data[dx])-1,cx,1,hit,0,15,15]) + '\n'
        output = IrysObj.cmapHead(len(data)) + output
        return output

class FilterSensitivity(): # determines sensitivity of repeatfilter function
    def __init__(self,diversity = 'order'):
        self.R_C = {'order':1e3, 'family':1e4,'genus':1e5,'species':1e6}[diversity] # range for repeat check (determines reading frame)
        self.C_D = {'order':3, 'family':3,'genus':2,'species':2}[diversity] # hits chromosome diffusion. marker for repeated hits
        self.H_D = {'order':5, 'family':4,'genus':4,'species':3}[diversity] # Repeated hits in same region of reference. higher value indicates a repeat.
        self.O = {'order':3, 'family':5,'genus':3,'species':1}[diversity] # in a scaffold, if hits on chromosome less than sensitivity.O (percentage), it will be omited
        self.F = {'order':5, 'family':5,'genus':5,'species':5}[diversity] # Value for hit filtering. <F continuous hits per scaffold will be deleted in refchr function.

def RepeatFilter(hit_dict):  # in: one scaffold each time, {chr:[(Pos_Que, Pos_Ref)(...)(...)],chr:[(Pos_Que, Pos_Ref)(...)(...)]...}
    Rpt_list, Qpt_list = {}, {}
    for chr_name in hit_dict: # step 1, detect and remove likely repeat elements
        for dx in hit_dict[chr_name]:
            Rpt_list[dx[0]//sensitivity.R_C] = Rpt_list[dx[0]//sensitivity.R_C] + [chr_name] if dx[0]//sensitivity.R_C in list(Rpt_list) else [chr_name] 
            Qpt_list[dx[1]//sensitivity.R_C] = Qpt_list[dx[1]//sensitivity.R_C] + [chr_name] if dx[1]//sensitivity.R_C in list(Qpt_list) else [chr_name]
    for dx in list(Rpt_list): # if hit can be found in 3 or more chromosome's, remove hit
        if len(set(Rpt_list[dx])) < sensitivity.C_D or len(Rpt_list[dx]) < sensitivity.H_D : del Rpt_list[dx]
    for dx in list(Qpt_list): # if hit goes to 4 different locations on the same chromosome, remove hit
        if len(Qpt_list[dx]) < sensitivity.H_D : del Qpt_list[dx]
    for chr_name in hit_dict: # remove everything in the repeat list present in original data
        hit_dict[chr_name] = [dx for dx in hit_dict[chr_name] if dx[0]//sensitivity.R_C not in list(Rpt_list)]
        hit_dict[chr_name] = [dx for dx in hit_dict[chr_name] if dx[1]//sensitivity.R_C not in list(Qpt_list)]
    for chr_name in list(hit_dict): # remove every empty object in dictionary 
        if hit_dict[chr_name] == [] : del hit_dict[chr_name]
    total_hit = sum(len(dx) for dx in hit_dict.values()) # step 2, rare non-repeat offtargets (few hits refer to a ref chromosome, remove these (it is due to random hits))
    for chr_name in filter(lambda dx: len(hit_dict[dx])*100/total_hit < sensitivity.O ,list(hit_dict)): del hit_dict[chr_name]
    return hit_dict

def HitGrouping(hit_list, noise_suppress = False): # hitlist = [(Pos_Que, Pos_Ref)(Pos_Que, Pos_Ref)...(Pos_Que, Pos_Ref)]
    ## fist hit in list evaluation:
    # check if hit is alone in hitlist (= no second hit available):
    if len(hit_list)<2: 
        # discard hit if noise_suppress is true (=insignificant)
        if noise_suppress == True:
            return []
        # Equal second hit to first hit: [ref_lonely_hit_start, ref_lonely_hit_end (=start), scaffold_lonely_hit_start, scaffold_lonely_hit_end (=start), count=1]
        else:
            return [hit_list[0][1],hit_list[0][1],hit_list[0][0],hit_list[0][0],'count=1']
    ## Determine where break must be introduced based on distance between two neighbouring hits:
    # sort scaffold hits based on their location on the scaffold:
    sorted_scaf_hits = sorted([hit[0] for hit in hit_list])
    # break distance = median of distance between sorted hits, multiplied with 100 => will be used as threshold for introducing break
    hit_index_list = range(1, len(sorted_scaf_hits))
    break_distance = sorted([sorted_scaf_hits[hit_index] - sorted_scaf_hits[hit_index-1] for hit_index in hit_index_list])[len(hit_index_list)//2-1] * 100

    ## Group hits based on break_distance evaluation:
    groups = []
    hit_count= 0
    # sort hits based on hit location on Ref:
    hit_sort = sorted(hit_list, key = lambda dx: int(dx[1])) 
    Q_start, R_start = hit_sort[0]
    ref_density_sum = 0
    for index in range(1,len(hit_list)-1):
        hit_count += 1 
        if abs(hit_sort[index][0] - hit_sort[index-1][0]) > break_distance or abs(hit_sort[index][1] - hit_sort[index-1][1]) > break_distance:
            R_end, Q_end = hit_sort[index-1][1], hit_sort[index-1][0]
            if hit_count > sensitivity.F:
                if noise_suppress == True and abs(R_start - R_end) / hit_count < 1000000:
                    groups.append([R_start, R_end, Q_start, Q_end, 'count='+str(hit_count)])
                elif noise_suppress == False:
                    groups.append([R_start, R_end, Q_start, Q_end, 'count='+str(hit_count)])
            Q_start, R_start = hit_sort[index+1]
            hit_count = 1

    ## When reaching end of hits in hit_list
    if noise_suppress == False or hit_count > sensitivity.F:
        # stl = second-to-last (hit)
        Q_end_stl,R_end_stl = hit_sort[-2][0], hit_sort[-2][1]
        Q_end_last, R_end_last = hit_sort[-1][0], hit_sort[-1][1]

        if noise_suppress == False:
            if abs(Q_end_last - Q_end_stl) > break_distance or abs(R_end_last - R_end_stl) > break_distance:
                groups += [[R_start, R_end_stl, Q_start, Q_end_stl,'count='+str(hit_count+1)]]
            else: groups += [[R_start, R_end_last, Q_start, Q_end_last,'count='+str(hit_count+1)]]
            
        else:
            if abs(R_start - R_end_last) / hit_count < 1000000:
                if abs(Q_end_last - Q_end_stl) > break_distance or abs(R_end_last - R_end_stl) > break_distance:
                    groups += [[R_start, R_end_stl, Q_start, Q_end_stl,'count='+str(hit_count+1)]]
                else: groups += [[R_start, R_end_last, Q_start, Q_end_last,'count='+str(hit_count+1)]]
    return groups # out: [[Ref_start,Ref_end,Que_start,Que_end,count=xxx][Ref_start,Ref_end,Que_start,Que_end,count=xxx]...], sorted on Ref

def ScaffoldAlignment(groups, inp, chrom):
    # section 1: sort scaffolds based on ref start of first section of list
    # groups: {scaffold:[[Ref_start,Ref_end,Que_start,Que_end,count=xxx], scaffold:[...]...}
    sorted_list, grouped_scaffolds, result_list, all_ref_start_ends = [], [], [], []
    for scaffold in groups.keys():
        sorted_list.append((scaffold, groups[scaffold][0][0]))
        all_ref_start_ends.append((scaffold, groups[scaffold][0][0], 'start_of_closest_con')) 
        all_ref_start_ends.append((scaffold, groups[scaffold][0][1], 'end_of_closest_con'))
    sorted_list.sort(key=operator.itemgetter(1))
    sorted_scaffolds = list(zip(*sorted_list))[0]
    current_scaffold = sorted_scaffolds[0]
    end_of_chr_flag, end_flag = False, False
    while end_of_chr_flag == False:
        ### run scaffold evaluation
        scaffold_evaluation = ComboCheck(inp, chrom, current_scaffold)
        
        ### calculate closest scaffold (and the side of the scaffold) to the current scaffold on the reference chromosome
        # remove current scaffold permanently out list of possibilities
        all_ref_start_ends = [i for i in all_ref_start_ends if i[0] != current_scaffold] 
        # end_flag checks on what side of the current_scaffold the next scaffold will be build on
        if end_flag == True: 
            current_pos = groups[sorted_scaffolds[sorted_scaffolds.index(current_scaffold)]][0][0] 
        else: current_pos = groups[sorted_scaffolds[sorted_scaffolds.index(current_scaffold)]][0][1]
        # calculate closest scaffold
        distance_list = []
        for con_name, ref_pos, side in all_ref_start_ends:
            distance_list.append((con_name, abs(current_pos - ref_pos), side))
        closest_distance = min(distance_list,key=operator.itemgetter(1))
        
        ### Update end_flag
        if closest_distance[2] == 'end_of_closest_con': end_flag = True
        else: end_flag = False
        
        ### Determine grouping of scaffolds
        grouped_scaffolds, result_list = Evaluation(scaffold_evaluation, current_scaffold, closest_distance, grouped_scaffolds, result_list)
        current_scaffold = closest_distance[0] # update scaffold to closest next one
        if len(all_ref_start_ends) == 2: 
            end_of_chr_flag = True # last remaining scaffold reached (two sides on one scaffold)
            scaffold_evaluation = ComboCheck(inp, chrom, current_scaffold)
            closest_distance = ('chr_end', 1, 'start_of_closest_con')
            grouped_scaffolds, result_list = Evaluation(scaffold_evaluation, current_scaffold, closest_distance, grouped_scaffolds, result_list)
            if grouped_scaffolds != []: result_list.append(grouped_scaffolds) 
    
    print('grouping: {}'.format(result_list))
    return result_list

def Evaluation(scaffold_evaluation, current_scaffold, closest_distance, grouped_scaffolds, result_list):
    # print(current_scaffold + ': ' + scaffold_evaluation)
    if scaffold_evaluation == 'right':
        if closest_distance[2] == 'start_of_closest_con':
            if grouped_scaffolds != []: result_list.append(grouped_scaffolds)
            grouped_scaffolds = []
            grouped_scaffolds.append(current_scaffold)
        else:
            grouped_scaffolds.append(current_scaffold)
            result_list.append(grouped_scaffolds)
            grouped_scaffolds = []
        
    elif scaffold_evaluation == 'left':
        if closest_distance[2] == 'end_of_closest_con':
            if grouped_scaffolds != []: result_list.append(grouped_scaffolds)
            grouped_scaffolds = []
            grouped_scaffolds.append(current_scaffold)
        else: 
            grouped_scaffolds.append(current_scaffold)
            result_list.append(grouped_scaffolds)
            grouped_scaffolds = []
        
    elif scaffold_evaluation == 'continue':
        grouped_scaffolds.append(current_scaffold)
        if closest_distance[1] > 25000000: 
            # print('break due to distance {}'.format(closest_distance))
            result_list.append(grouped_scaffolds)
            grouped_scaffolds = []

    elif scaffold_evaluation == 'stop':
        if grouped_scaffolds != []: result_list.append(grouped_scaffolds)
        grouped_scaffolds = []
        result_list.append([current_scaffold])
    
    return grouped_scaffolds, result_list

def ComboCheck(inp, current_chr, scaffold_to_check):
    # input {scaffold:{chrX:[[Ref_start,Ref_end,Que_start,Que_end,count=xxx],...],chrX: [...]}}, with list sorted on Ref
    for scaffold in inp.keys():
        if scaffold == scaffold_to_check:
            # Check if scaffold is composed out of more than one ref chromosome section
            chr_count = len(inp[scaffold].keys())
            if chr_count == 1: evaluation = 'continue'
            else:
                ref_min, ref_max, ref_min_cur, ref_max_cur = 0, 0, 0, 0
                for chrom in inp[scaffold].keys():
                    block = inp[scaffold][chrom] # list of values belonging to certain chrom
                    for section in block: # for each sublist in list
                        if chrom == current_chr: # check chromosome 
                            # avoid possible issues
                            if ref_min == 0 or ref_max == 0:
                                ref_min_cur, ref_min = section[2], section[2]
                                ref_max_cur, ref_max = section[3], section[3]
                            if ref_min_cur == 0 or ref_max_cur == 0: ref_min_cur, ref_max_cur = section[2],section[3]
                            
                            if section[2] < ref_min_cur:
                                ref_min_cur = section[2]
                                if ref_min > ref_min_cur: ref_min = ref_min_cur
                            if section[3] > ref_max_cur:
                                ref_max_cur = section[3]
                                if ref_max < ref_max_cur: ref_max = ref_max_cur
                        
                        else:
                            if ref_min == 0 or ref_max == 0: ref_min, ref_max = section[2], section[3]
                            if ref_min > section[2]: ref_min = section[2]
                            if ref_max < section[3]: ref_max = section[3]
                # run check
                flagL, flagR, count = True, True, 0
                if ref_min != ref_min_cur:
                    count += 1
                    flagL = False
                if ref_max != ref_max_cur:
                    count += 1
                    flagR = False

                if count == 2:
                    evaluation = 'stop' # current_chr is isolated
                elif count == 1: 
                    if flagL == False: evaluation = 'right' # scaffold can be attached to other scaffolds on the right
                    if flagR == False: evaluation = 'left' # scaffold can be attached to other scaffolds on the left
                else:
                    evaluation = 'continue' # both sides of the combination consist out of the current_chr
    return evaluation

def ContigGrouping(hit_dict,scaffold_name):
    group_list = {}
    for dx in hit_dict[scaffold_name].keys(): group_list[dx] = HitGrouping(hit_dict[scaffold_name][dx], noise_suppress = True)
    return group_list

def QueryFaiOpen(Call_from):
    if not args.faindex: 
        print('Fasta index file is needed for '+Call_from+' function.')
        exit(1)
    with open(args.faindex,'r') as FAI: 
        return {dx.split('\t')[0]:int(dx.split('\t')[1]) for dx in FAI.readlines()}

def output_files(filename, height, box_list, line_list, text_list, pdf):
    if filename[-4:] == ".svg" or filename[-4:] == ".pdf": filename = filename[:-4]

    ## create svg file
    with open(filename + '.svg','w') as SVG_out: SVG_out.write(SVG.head(2000, height) + '\n  '.join(box_list)+ '\n  '.join(line_list) + '\n  '.join(text_list) + '</svg>\n')
    if pdf == True:
        ## Wand module converts svg to png, than img2pdf is used to convert png to pdf
        with Image(filename=filename + '.svg') as img:
            img.alpha_channel = False
            img.format = 'png'
            img.save(filename='temp.png')
        os.remove(filename + '.svg')
        with open(filename + '.pdf',"wb") as pdf: pdf.write(img2pdf.convert('temp.png'))
        os.remove('temp.png')

############################################## argparse module ################################################

cwd = os.getcwd()
dir_name = cwd.split('/')
dir_name = dir_name[-1]

parser = argparse.ArgumentParser(prog='SynChroBuild', description='Build synteny maps from BLAST results (Synteny Builder) and output the synteny data into AGP, fasta or summary files. Uses python 3.4+')
parser.add_argument('blastn', help = 'Input file name of BLASTN result (for synteny map) or fasta file (for AGP, fasta and summary files)')
parser.add_argument('faindex', help = 'Input fasta index file (for synteny map) or synteny file (for AGP, fasta and summary files)')

# Blast file input
parser_syn2 = parser.add_argument_group('Synteny builder - main functions:')
parser_syn2.add_argument('-s','--scaf', dest='scaffold_name', action ='store', default = '', metavar = 'SCAF', help = 'Generate SVG file with detailed blastn hits of a specified scaffold')
parser_syn2.add_argument('-c','--chr', dest='ref_chr', action ='store', nargs = '?', const = 'All', default = False, metavar='CHROM', help = 'Build synteny on chosen chromosome (default = all chromosomes).')
parser_syn2.add_argument('-a','--align', dest='scaf_aligner', action ='store_true', default = False, help = 'Align contigious scaffolds into groups in the SVG output file.')

parser_syn3 = parser.add_argument_group('Synteny builder - secundary functions:')
parser_syn3.add_argument('-x','--xmap', dest='Xmap_function', action ='store_true', help = 'Output list as Irysview .cmap and .xmap format.')
parser_syn3.add_argument('-g','--g_band', dest='G_scaffold', action ='store', default = '', metavar = 'SCAF', help='Generating G-banding simulation of scaffold based on synteny. Output .SVG format.') 
parser_syn3.add_argument('-b','--blast_sum', dest='blastn_summary', action='store_true', help='Give statistics about the BLASTN result and list all significant scaffolds')
# create text file = density of blast (hit ratio), simularity between both species, average

parser_syn1 = parser.add_argument_group('Synteny builder - settings:')
parser_syn1.add_argument('-r','--ref',dest='Ref_genome',action = 'store',default = 'MM10', metavar='FAIDEX', help='Information of Reference genome loaded from a fasta index file, default = GRCm38/MM10.')
parser_syn1.add_argument('-d','--div',dest='diversity',action ='store',default = 'order', metavar='VALUE', help='Set diversity value between Query and Reference(order/family/genus/species)')
parser_syn1.add_argument('-m','--offset', dest='offset', action ='store', default = '0', metavar='BP', help = 'Shift the scaffold to right (+) or left (-) compared to the reference chromosome. Units in basepair, K,M,G also accepted. Example: -m"-1000" = -m"-1K"')
parser_syn1.add_argument('-e','--exp', dest='exp_value', action ='store', default = '200', metavar='VALUE', help = 'Set threshold for expectation (E) value of BLAST results, e.g. -e 80 indicate E value < 1e-80.')
parser_syn1.add_argument('-o', dest='output', action = 'store', metavar='NAME', nargs='?', const = dir_name, default = dir_name, help = 'Assign a custom name to the output file(s).')
parser_syn1.add_argument('-pdf', dest='pdf_flag', action ='store_true', default = False, help='Output file in pdf format instead of svg format.')
args = parser.parse_args()

############################################## script ################################################

### Check input parameters
## Check for excessive arguments
if sum([args.blastn_summary, args.ref_chr!= False, args.scaffold_name != '', args.Xmap_function], args.G_scaffold != '') != 1:
    parser.print_help()
    exit(1)

if args.blastn[-4:] == '.fai': # input has been switched by accident
    faindex = args.blastn
    args.blastn = args.faindex
    args.faindex = faindex

if args.Ref_genome != 'MM10':
    genome = args.Ref_genome
else: genome = 'MM10'

## Check offset section (-m): Moves scaffold (requires -f option)
input_offset = args.offset
# if input is a positive number
if input_offset.isdigit(): offset = int(input_offset)/1e5
# if input is a negative number 
elif input_offset[0] == '-' and input_offset[1:].isdigit(): offset = int(input_offset)/1e5
# if input uses K, M or G and is positive
elif input_offset[-1] in ['K','M','G'] and input_offset[:-1].isdigit(): offset = int(input_offset[:-1])*{'K':1e-2,'M':10,'G':1e4, 'k':1e-2,'m':10,'g':1e4}[input_offset[-1]]
# if input uses K, M or G and is negative 
elif input_offset[-1] in ['K','M','G'] and input_offset[0] == '-' and input_offset[1:-1].isdigit(): offset = int(input_offset[1:-1])*{'K':1e-2,'M':10,'G':1e4, 'k':1e-2,'m':10,'g':1e4}[input_offset[-1]]
else:
    print('--offset argument: must be a number, or end with K,M,G.')
    exit(1)

### Create hit_list dictionary fron blast input
# structure of hit_list: hit_list = {scaffold:{chr:[(scaffold_position, chr_position)(...)(...)]} } => chr_pos is centre of blast hit location on both ref and scaffold
hit_list = {} 

# set repeats filter sensitivity + set threshold E-value for BLAST results
sensitivity = FilterSensitivity(parser.parse_args().diversity) 
expection = int(parser.parse_args().exp_value)

hit_sort, data = [], []
with open(args.blastn,'r') as BLAST: 
    while True:
        data = BLAST.readline().split('\t')
        if data == [''] : break
        scaffold_name, chr_name = data[0], data[1]
        scaffold_start, match_start, match_end, scaffold_end, E_value = int(data[6]), int(data[7]), int(data[8]), int(data[9]), data[10]
        if E_value == '0.0' or (E_value.split('-')[0] == '1e' and int(E_value.split('-')[1]) == expection) or int(E_value.split('-')[1]) > expection:
            # if scaffold already exists in hit_list
            if scaffold_name in hit_list.keys(): 
                if chr_name.lower() in hit_list[scaffold_name].keys(): 
                    hit_list[scaffold_name][chr_name.lower()] += [tuple([(scaffold_start + match_start)//2,(match_end + scaffold_end)//2])] 
                else: hit_list[scaffold_name][chr_name.lower()] = [tuple([(scaffold_start + match_start)//2,(match_end + scaffold_end)//2])]
            # if scaffold is not yet present in hit_list
            else:
                hit_list[scaffold_name] = {chr_name.lower():[tuple([(scaffold_start + match_start)//2,(match_end + scaffold_end)//2])]}
        # update list of chromosomes if not yet added for all chromosomes

### Filter blast result: Remove repeats
for scaffold_name in hit_list.keys(): hit_list[scaffold_name] = RepeatFilter(hit_list[scaffold_name])

### Create summary of BLASTN result
if args.blastn_summary:
    # BLASTN hit count for each scaffold/chr
    simularity_sum, hit_length_sum, gap_length_sum, average_count = 0, 0, 0, 0
    chr_hit_count, scaffold_count = 0, 0
    prev_scaffold_name, prev_chr_name = '', ''
    chr_hit_dict, scaffold_dict = {}, {}
    BLAST = open(args.blastn, 'r')
    for line in BLAST.readlines():
        line = line.split('\t')
        scaffold_name, chr_name, simularity, hit_length, gap_length = line[0], line[1], line[2], line[3], line[4]
        if scaffold_name in hit_list.keys(): # Filter scaffolds
            # statistics
            simularity_sum += float(simularity)
            hit_length_sum += int(hit_length)
            gap_length_sum += int(gap_length)
            average_count += 1
            # scaffold data
            if prev_scaffold_name == '': prev_scaffold_name = scaffold_name
            if scaffold_name != prev_scaffold_name:
                chr_hit_dict[prev_chr_name.lower()] = chr_hit_count
                scaffold_dict[prev_scaffold_name] = chr_hit_dict
                chr_hit_dict = {}
                scaffold_count += 1
                chr_hit_count = 0
            if chr_name != prev_chr_name:
                chr_hit_dict[prev_chr_name.lower()] = chr_hit_count
                chr_hit_count = 0
            chr_hit_count += 1
            prev_scaffold_name, prev_chr_name = scaffold_name, chr_name
    
    # counting chromosomes hit
    scaffold_length = QueryFaiOpen('ListAll')
    result_list, total_length = [], 0
    for scaffold_name in hit_list.keys():
        if scaffold_name in scaffold_dict.keys():
            result_list.append(str(scaffold_name + '\t' + str(scaffold_length[scaffold_name]) + '\t' + str(sum(scaffold_dict[scaffold_name].values())) + '\n'))
            total_length += scaffold_length[scaffold_name]
            for chr_name in hit_list[scaffold_name].keys():
                hitgroup = HitGrouping(hit_list[scaffold_name][chr_name])
                for value in hitgroup: # [[1,2,3,4,5],[1,2,3,4,5]] or [1,2,3,4,5]
                    if type(value) == list:
                        hitgroup = value # [1,2,3,4,5]
                        result_list.append(chr_name + '\t' + str(scaffold_dict[scaffold_name][chr_name]) + '\t' + str(hitgroup[0]) + '\t' + str(hitgroup[1]) + '\t' + str(hitgroup[2]) + '\t' + str(hitgroup[3]) + '\n')
                    else: 
                        result_list.append(chr_name + '\t' + str(scaffold_dict[scaffold_name][chr_name]) + '\t' + str(hitgroup[0]) + '\t' + str(hitgroup[1]) + '\t' + str(hitgroup[2]) + '\t' + str(hitgroup[3]) + '\n')
                        break

    # statistics
    simularity_file = str(round(simularity_sum/average_count,2)) + '%'
    hit_length_file = str(round(hit_length_sum/average_count)) + ' bp'
    gap_length_file = str(round(gap_length_sum/average_count)) + ' bp'
    statistics = '### BLASTN statistics ###\nAmount of scaffolds: ' + str(scaffold_count) + '\n' + 'Combined scaffold length: ' + str(total_length) + ' bp\n' + 'Average simularity: ' + simularity_file + '\n' + 'Average hit length: ' + hit_length_file + '\n' + 'Average gap length: ' + gap_length_file + '\n\n'
    
    output_file = open(dir_name + '_BLASTN_summary.txt','w')
    output_file.write(statistics)
    output_file.write('### list of filtered scaffolds ###\n# scaffold   scaffold length   blast hits\n# chr   hits/chr   [Ref_chr_start, Ref_chr_end, Target_chr_start, Target_chr_end]\n# ...')
    print(statistics)
    for line in result_list:
        output_file.write(line)
    output_file.close()
    exit(0)

### Generate a SVG file for one specific scaffold
if args.scaffold_name != '':

    reverse = False
    ## check if an inverted scaffold orientation is needed:
    all_scaf_hits, all_chr_hits = [], []
    for chrom in hit_list[args.scaffold_name].keys():
        for hit in hit_list[args.scaffold_name][chrom]:
            all_scaf_hits.append(hit[0])
            all_chr_hits.append(hit[1])
    scaf_median = statistics.median(all_scaf_hits)
    stored_index = []
    ## filter hit outliers
    for hit in all_scaf_hits:
        if hit < scaf_median + 3000000 and hit > scaf_median - 3000000: stored_index.append(all_scaf_hits.index(hit))
    all_scaf_hits = [all_scaf_hits[i] for i in stored_index]
    all_chr_hits = [all_chr_hits[i] for i in stored_index]
    min_chr_index = all_scaf_hits.index(min(all_scaf_hits))
    max_chr_index = all_scaf_hits.index(max(all_scaf_hits))
    if all_chr_hits[max_chr_index] < all_chr_hits[min_chr_index]: reverse = True

    ## Generate reference chromosome svg data
    scaffold_length = QueryFaiOpen('Figure')
    SVG = SvgObj()
    line_list = []
    chr_median = statistics.median(all_chr_hits)
    if chr_median > 175000000: chr_median = 175000000
    # defines scaffold name printed in svg file
    box_print, text_print = SVG.ChrBox(round(chr_median/1e5 + int(offset)),300, args.scaffold_name, genome, False, optionP = False, chr_detail = [scaffold_length[args.scaffold_name],[]])
    # lists will be printed in the end during svg creation (box_print contains svg path object data, text_print contains text data)
    box_list, text_list = [box_print], [text_print]
    total_hit = sum(len(dx) for dx in hit_list[args.scaffold_name].values())
    hit_sort = sorted(zip([dx for dx in hit_list[args.scaffold_name].keys()],[len(hit_list[args.scaffold_name][dx]) for dx in hit_list[args.scaffold_name].keys()]), key = lambda dx: dx[1], reverse = True)
    
    ## Generate scaffold svg data
    count = 0.0
    for cx,chr_name in enumerate(hit_sort):
        if cx == 0: cx = -1
        # defines chr name printed in svg file
        box_print, text_print = SVG.ChrBox(0,300+cx*200,chr_name[0], genome, True, optionP = False, chr_detail=[])
        box_list.append(box_print)
        text_list.append(text_print)
        for rec in hit_list[args.scaffold_name][chr_name[0]]:
            # Slightly different colours compared to the colours of reference chromosomes
            hit_col_dict = {'chr1':'#F06565','chr2':'#F77D2C','chr3':'#EEBE00','chr4':'#A3CB00','chr5':'#56D800','chr6':'#00D0A6','chr7':'#00005C','chr8':'#2E78E6','chr9':'#8F4EF2','chr10':'#EFA6E0','chr11':'#CE0505','chr12':'#652800','chr13': '#846C0B','chr14':'#A28512','chr15': '#428714','chr16': '#082B08','chr17': '#004CBD','chr18': '#132D54','chr19': '#B20F91', 'chr20':'#C67D00', 'chr21':'#150804','chr22':'#C3C354','chr23':'#3F3F1C','chr24':'#1C4A3E','chr25':'#2FD4B8','chr26':'#D11ABF','chr27':'#4A0743','chr28':'#803578','chr29':'#958EC0','chr30':'#5B5E01','chr31':'#A5B779','chr32':'#BF7B48','chr33':'#00DB66','chr34':'#DEE43A','chr35':'#4D0000','chr36':'#9476C0','chr37':'#805F25','chr38':'#D15095','chrx': '#2C016E','other': '#7B797D'}
            if cx == -1: # if section of scaffold contains the same chrom as the ref chrom
                if reverse == True: line_list.append('<path d= "M ' + str(int(rec[1])/1e5) + ' ' + str(300 + cx*200) +' L ' + str(int(rec[1])/1e5) + ' ' + str(400+ cx*200) + ' L ' + str(int(int(scaffold_length[args.scaffold_name])/1e5 - int(rec[0])/1e5 + chr_median/1e5 + offset)) + ' 300 L ' + str(int(int(scaffold_length[args.scaffold_name])/1e5 - int(rec[0])/1e5 + chr_median/1e5 + offset)) + ' 400" fill = "none" stroke = "' + hit_col_dict[chr_name[0]] + '" stroke-width="1" />')
                else: line_list.append('<path d= "M ' + str(int(rec[1])/1e5) + ' ' + str(300 + cx*200) +' L ' + str(int(rec[1])/1e5) + ' ' + str(400+ cx*200) + ' L ' + str(int(int(rec[0])/1e5 + chr_median/1e5 + offset)) + ' 300 L ' + str(int(int(rec[0])/1e5 + chr_median/1e5 + offset)) + ' 400" fill = "none" stroke = "' + hit_col_dict[chr_name[0]] + '" stroke-width="1" />')
            else: # if section of scaffold contains a different chrom than the ref chrom
                if reverse == True: line_list.append('<path d= "M ' + str(int(int(scaffold_length[args.scaffold_name])/1e5 - int(rec[0])/1e5 + chr_median/1e5 + offset)) + ' 300 L ' + str(int(int(scaffold_length[args.scaffold_name])/1e5 - int(rec[0])/1e5 + chr_median/1e5 + offset)) + ' 400 L ' + str(int(rec[1])/1e5) + ' ' + str(300 + cx*200) +' L ' + str(int(rec[1])/1e5) + ' ' + str(400+ cx*200) +'" fill = "none" stroke = "' + hit_col_dict[chr_name[0]] + '" stroke-width="1" />')
                else: line_list.append('<path d= "M ' + str(int(int(rec[0])/1e5 + chr_median/1e5 + offset)) + ' 300 L ' + str(int(int(rec[0])/1e5) + chr_median/1e5 + offset) + ' 400 L ' + str(int(rec[1])/1e5) + ' ' + str(300 + cx*200) +' L ' + str(int(rec[1])/1e5) + ' ' + str(400+ cx*200) +'" fill = "none" stroke = "'+ hit_col_dict[chr_name[0]] + '" stroke-width="1" />')
        count = count + len(hit_list[args.scaffold_name][chr_name[0]])
        if count / total_hit > 0.95: 
            height = 800 + cx*200
            break
    line_list = list(set(line_list))
    
    ## output
    if args.output != '':
        filename = args.output + '_' + args.scaffold_name
        output_files(filename, height, box_list, line_list, text_list, args.pdf_flag)
    else: print(SVG.head(2000, len(height)*100+700) + '\n  ' + '\n  '.join(list(box_list)) + '\n  '.join(line_list) + '\n  '.join(text_list) + '</svg>\n')
    exit(0)

### Generate a svg file containing all scaffolds of a reference chromosome.
if args.ref_chr != False:
    optionP = False # optionP is a flag used to indicate the usage of the scaffold align option, also used as counter if True.
    chr_list = []
    if args.ref_chr != 'All': 
        chr_name = args.ref_chr.lower()
        chr_list.append(chr_name)
    else:
        ## create list of all chromosomes
        chr_list = []
        for scaffold_name in hit_list.keys():
            for chromosome in hit_list[scaffold_name].keys():
                in_list = chromosome.lower() in chr_list
                if in_list == False: chr_list.append(chromosome)
            if 'chry' in chr_list: 
                chr_list[chr_list.index('chry')] = 'other'

    ## hitgrouping (group hits per scaffold)
    # INPUT HitGrouping: hit_list: {scaffold:{chr:[(Que, Ref)(Que, Ref)(...)(...)],chr:[(...)(...)(...)]}}
    # OUTPUT HitGrouping: groups: {ctgname:[[Ref_start,Ref_end,Que_start,Que_end,count=xxx],[Ref_start,Ref_end,Que_start,Que_end,count=xxx]]...}
    for chr_name in chr_list:
        scaffold_length = QueryFaiOpen('Synteny')
        SVG = SvgObj()
        # Determines ref chromosome and name
        box_print, text_print = SVG.ChrBox(0, 300, chr_name, genome, False, optionP)
        box_list, text_list = [box_print], [text_print]
        line_list, boundary_list, groups = [], [], {} 
        
        for scaffold_name in hit_list.keys():
            if chr_name in hit_list[scaffold_name].keys():
                groups[scaffold_name] = HitGrouping(hit_list[scaffold_name][chr_name], noise_suppress = True)
                for dx in hit_list[scaffold_name].keys():
                    if len(hit_list[scaffold_name][dx]) > sensitivity.O and dx != chr_name and scaffold_name not in boundary_list:
                        boundary_list.append(scaffold_name)
        groups = dict((dx, ex) for dx, ex in groups.items() if ex != [])
        
        ## Section specific for obtaining pseudochromosome data
        if args.scaf_aligner != False:
            input_dictionary = {}
            for dx in groups:
                chr_specific_scaffolds = {}
                for chr_ in hit_list[dx].keys():
                    block_list = []
                    chr_block = HitGrouping(hit_list[dx][chr_], noise_suppress = True)
                    if chr_block != []:
                        for block in chr_block:
                            block_list.append(block)
                        #dict1 = {chr_:block_list}
                        chr_specific_scaffolds.update({chr_:block_list})
                input_dictionary.update({dx:chr_specific_scaffolds})
            result_list = ScaffoldAlignment(groups, input_dictionary, args.ref_chr)
            
            # create index_list that gives scaffolds in the same group the same id
            indexnumber, indexed_list = 1, []
            while indexnumber <= len(result_list):
                for scaffold in result_list[indexnumber - 1]:
                    indexed_list.append([indexnumber, scaffold])
                indexnumber += 1
            index = indexed_list

        else: index = enumerate(groups)

        ## Set position parameters of each object in svg file and create svg object
        # Set position parameters:
        previous_cx = 0
        added_length = 0
        counter = 1
        for cx, dx in index: # create index of remaining scaffolds
            baseline, bottomline = 400, 0
            if args.scaf_aligner != False: 
                optionP = counter
                if cx != previous_cx: 
                    # added_length = 0
                    chr_hit_list = []
                    for hit in hit_list[dx][chr_name.lower()]:
                        chr_hit_list.append(hit[1])
                    chr_median = statistics.median(chr_hit_list) / 1e5
                    if chr_median > 160000000: chr_median = 160000000
                    added_length = chr_median - 100
                xco = round(added_length)
                if cx % 3 == 0: yco = 600
                elif cx % 4 == 0: yco, baseline, bottomline = 10, 300, 100
                elif cx % 2 == 0: yco, baseline, bottomline = 120, 300, 100
                else: yco = 800
                
                # update loop values
                previous_cx = cx
                added_length += scaffold_length[dx]/1e5 + 5
                counter += 1
                
            else:
                xco = int((groups[dx][0][0]-min(groups[dx][0][2:4]))/1e5)
                yco = 500 + cx*100

            # Create svg object
            # determines scaffold box and name
            box_print, text_print = SVG.ChrBox(xco, yco, str(dx), genome, False, optionP, chr_detail = [scaffold_length[dx],[]])
            box_list.append(box_print)
            text_list.append(text_print)

            previous_chr = ''
            if dx in boundary_list:
                for chr_ in hit_list[dx].keys():
                    chr_block = HitGrouping(hit_list[dx][chr_], noise_suppress = True)
                    if chr_block != []:
                        for _,_,start,end,_ in chr_block:
                            if previous_chr != chr_:
                                # Create chr blocks for each scaffold
                                box_print, text_print = SVG.ChrBox(xco + (min(start,end)/1e5),yco, chr_, genome, False, optionP, chr_detail=[abs(end-start),[]])
                                box_list.append(box_print) 
                                text_list.append(text_print)
                            else:
                                box_print,_ = SVG.ChrBox(xco + (min(start,end)/1e5),yco, chr_, genome, False, optionP, chr_detail=[abs(end-start),[]])
                                box_list.append(box_print)

                            # update loop values
                            previous_chr = chr_

            # creates the connecting figure between the ref chrom and the scaffold 
            for ex in range(len(groups[dx])):
                corner1 = '"M'+str(int(groups[dx][ex][0]/1e5)) + ' ' + str(baseline)
                corner2 = ' L'+str(int(groups[dx][ex][1]/1e5)) + ' ' + str(baseline)
                corner3 = ' L'+str(xco + (int(groups[dx][ex][3])/1e5)) + ' ' + str(yco + bottomline)
                corner4 = ' L'+str(xco + (int(groups[dx][ex][2])/1e5)) + ' ' + str(yco + bottomline) + ' Z"'
                line_list.append('<path d='+corner1+corner2+corner3+corner4+' fill="black" stroke="none" fill-opacity="0.4" />')
    
            ## output
            if args.output != '':
                filename = args.output
                if filename[-4:] == ".svg" or filename[-4:] == ".pdf":
                    filename = filename[:-4]
                filename = filename + '_' + chr_name
                output_files(filename, len(groups)*100+700, box_list, line_list, text_list, args.pdf_flag)
            else: print(SVG.head(2000, len(groups)*100+700) + '\n  ' + '\n  '.join(list(box_list)) + '\n  '.join(line_list) + '\n  '.join(text_list) + '</svg>\n')

### G-banding section
if args.G_scaffold != '':
    scaffold_length = QueryFaiOpen('G-banding')
    group_list = {}
    Band_SVG = ChrBand()
    SVG = SvgObj()
    for chr_ in hit_list[args.G_scaffold].keys():
        # group_list: {chr:[[Ref_start,Ref_end,Que_start,Que_end,count=xxx][Ref_start,Ref_end,Que_start,Que_end,count=xxx]],chr:[...]}
        group_list[chr_] = HitGrouping(hit_list[args.G_scaffold][chr_], noise_suppress=True)
    
    if args.output != '':
        with open(args.output,'w') as SVG_out:
            SVG_out.write(SVG.head(int(scaffold_length[args.G_scaffold]/1e6/7*100), 500) + '\n  ' + Band_SVG.BandSVG(scaffold_length[args.G_scaffold],7,group_list,x=0,y=0,height=100)+ '</svg>\n')
    else:
        print(SVG.head(int(scaffold_length[args.G_scaffold]/1e6/7*100), 500) + '\n' + Band_SVG.BandSVG(scaffold_length[args.G_scaffold],7,group_list,x=0,y=0,height=100)+ '</svg>\n')
