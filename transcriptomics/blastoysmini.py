#!/usr/bin/env python
#Written by ZNR; 2014/09/07
import argparse
import curses
import copy
import os
import json
import subprocess
import sys
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def transrender(screen,transstrand,stranddex):
	#Reset cursor and view position
	basedex = 0
	readdex = 0
	xview = 0
	yview = 0
	#Fix graphics boundaries
	leftbound = 1
	upperbound = 6
	lowerbound=curses.LINES-1
	rightbound=curses.COLS-1
	#Misc vars
	curstrand = transstrand[stranddex]
	sellist = []
	while 1:
		#Calculate cursor locations
		basecursor = basedex-xview
		readcursor = readdex-yview
		numreads = len(curstrand.reads)
		curstrand.kend = curstrand.kbegin + curstrand.kwidth
		# Refresh screen and labels
		screen.clear()
		screen.border(0)
		screen.addstr(1, 1, "Strand ID:"+str(stranddex)+"of"+str(len(transstrand)-1))
		screen.addstr(1, 20, "Min depth:"+str(curstrand.mindepth))
		screen.addstr(1, 35, "kmer width:"+str(curstrand.kwidth))
		screen.addstr(2, 1, "q to quit, wasd to move view (shift to go faster), op to switch strand, arrows to set kmer bounds, -=[] to move selector (shift to go faster)")
		screen.addstr(3, 1, "r to generate reference,  m to switch min depth, n to copy to new transcript, x to remove read, shift-Z to delete transcript")
		# Blit
		for xblit in range(leftbound,rightbound):
			xindex = xblit-leftbound+xview
			#individual strand status

			#transcript
			if xindex == basedex: rcolor = 10
			elif xindex in range(curstrand.kbegin,curstrand.kend): rcolor = 5
			else: rcolor = 0
			basechar = str(curstrand.transcript[xindex:xindex+1])
			yblit = 5
			if 'A' in basechar: screen.addstr(yblit,xblit,'A',curses.color_pair(rcolor+1))
			elif 'C' in basechar: screen.addstr(yblit,xblit,'C',curses.color_pair(rcolor+2))
			elif 'G' in basechar: screen.addstr(yblit,xblit,'G',curses.color_pair(rcolor+3))
			elif 'T' in basechar: screen.addstr(yblit,xblit,'T',curses.color_pair(rcolor+4))
			else: screen.addstr(yblit,xblit,'n',curses.color_pair(rcolor))
			#reads
			for yblit in range(upperbound,lowerbound):
				yindex = yblit-upperbound+yview
				if yindex < numreads:
					yread = curstrand.reads[yindex]
					yalign = curstrand.aligns[yindex]
					yend = curstrand.ends[yindex]
					#Row Highlighting
					if yindex == readdex or xindex == basedex: rcolor = 10
					elif yindex in sellist: rcolor = 5
					else: rcolor = 0
					#Blit bits
					if xindex<yalign: screen.addstr(yblit,xblit,'.',curses.color_pair(rcolor))
					elif xindex>=yend: screen.addstr(yblit,xblit,'.',curses.color_pair(rcolor))
					elif yalign<=xindex<yend:
						basechar = str(yread[xindex-yalign:xindex-yalign+1])
						if 'A' in basechar: screen.addstr(yblit,xblit,'A',curses.color_pair(rcolor+1))
						elif 'C' in basechar: screen.addstr(yblit,xblit,'C',curses.color_pair(rcolor+2))
						elif 'G' in basechar: screen.addstr(yblit,xblit,'G',curses.color_pair(rcolor+3))
						elif 'T' in basechar: screen.addstr(yblit,xblit,'T',curses.color_pair(rcolor+4))
						else: screen.addstr(yblit,xblit,'X',curses.color_pair(rcolor))
		screen.refresh()
		# screen.nodelay(1)
		# Get and interpret keyboard input
		y = screen.getch()
		if y == ord('q'):	# Quitting options
			stranddex = -1
			break
		elif y == ord('o') and stranddex != 0:
			stranddex += -1
			break
		elif y == ord('p') and stranddex+1 < len(transstrand):
			stranddex += 1
			break
		elif y == ord('-') and basedex > 0: basedex +=-1	# Navigation options
		elif y == ord('_') and basedex > 5: basedex +=-5
		elif y == ord('=') and basedex < len(curstrand.transcript): basedex +=1
		elif y == ord('+') and basedex < len(curstrand.transcript)-5: basedex +=5
		elif y == ord('[') and readdex > 0: readdex +=-1
		elif y == ord('{') and readdex > 5: readdex +=-5
		elif y == ord(']') and readdex+1 < numreads: readdex +=1
		elif y == ord('}') and readdex+1 < numreads-5: readdex +=5
		elif y == ord('w') and yview>0: yview +=-1
		elif y == ord('W') and yview>5: yview +=-5
		elif y == ord('a') and xview>0: xview +=-1
		elif y == ord('A') and xview>5: xview +=-5
		elif y == ord('s') and yview<(numreads): yview +=1
		elif y == ord('S') and yview<(numreads)-5: yview +=5
		elif y == ord('d') and xview+(rightbound-leftbound)<len(curstrand.transcript): xview +=1
		elif y == ord('D') and xview+(rightbound-leftbound)<len(curstrand.transcript)-5: xview +=5
		elif y == curses.KEY_LEFT and curstrand.kbegin>0: curstrand.kbegin +=-1
		elif y == curses.KEY_RIGHT and curstrand.kwidth+curstrand.kbegin<len(curstrand.transcript): curstrand.kbegin  +=1
		elif y == curses.KEY_UP and curstrand.kwidth>0: curstrand.kwidth +=-1
		elif y == curses.KEY_DOWN and curstrand.kwidth+curstrand.kbegin<len(curstrand.transcript): curstrand.kwidth +=1
		elif y == ord(' '):	# Misc Options
			if readdex in sellist: sellist.remove(readdex)
			else: sellist.append(readdex)
		elif y == ord('n'):
			transstrand.append(copy.deepcopy(curstrand))
			transstrand[-1].name = str(len(transstrand)-1)+"_"+curstrand.name
		elif y == ord('Z'):
			transstrand.pop(stranddex)
			stranddex += -1
			break
		elif y == ord('x'):
			curstrand.reads.pop(readdex)
			curstrand.ids.pop(readdex)
			curstrand.aligns.pop(readdex)
			curstrand.ends.pop(readdex)
			for n in sellist:
				if n == readdex: sellist.remove(n)
				elif n > readdex:
					sellist.remove(n)
					sellist.append(n-1)
			readdex +=-1
			if readdex < 0: readdex=0
		elif y == ord('r'):
			newtranscript = ''
			for n in range(0,len(curstrand.transcript)):
				ctA=0
				ctC=0
				ctG=0
				ctT=0
				for o in range(0,numreads):
					yread = curstrand.reads[o]
					yalign = curstrand.aligns[o]
					yend = curstrand.ends[o]
					if yalign<=n<yend:
						basechar = str(yread[n-yalign:n-yalign+1])
						if basechar in 'A': ctA +=1
						elif basechar in 'C': ctC +=1
						elif basechar in 'G': ctG +=1
						elif basechar in 'T': ctT +=1
				maxN = max([ctA,ctC,ctG,ctT])
				if maxN == ctA and max([ctC-maxN,ctG-maxN,ctT-maxN]) <= 0-curstrand.mindepth: newtranscript += 'A'
				elif maxN == ctC and max([ctA-maxN,ctG-maxN,ctT-maxN]) <= 0-curstrand.mindepth: newtranscript += 'C'
				elif maxN == ctG and max([ctA-maxN,ctC-maxN,ctT-maxN]) <= 0-curstrand.mindepth: newtranscript += 'G'
				elif maxN == ctT and max([ctA-maxN,ctC-maxN,ctG-maxN]) <= 0-curstrand.mindepth: newtranscript += 'T'
				else: newtranscript += 'n'
			curstrand.transcript=Seq(newtranscript,generic_dna)
		elif y == ord('m'):
			depths = [1,2,3,4,5,10]
			curstrand.mindepth = depths[depths.index(curstrand.mindepth)-1]
	return stranddex

class BLASToys:
	def transextend(self, filepath, temppath, revcom):
		searchseed = self.transcript[self.kbegin:self.kend]
		aln1 = self.transcript.find(searchseed)
		txfound = 0
		# must have single line fastas to work - no n60 splitting; otherwise, can't use grep for max speed
		with open(temppath + '/temp','w') as tempfile:
			if revcom == True:
				grep = subprocess.Popen(['grep', '-B1', revcom(searchseed), filepath],stdout=tempfile)
			else:
				grep = subprocess.Popen(['grep', '-B1', searchseed, filepath],stdout=tempfile)
			grep.wait()
			tempfile.flush()
		with open(temppath + '/temp','r') as tempfile:
			line1=''
			line2=''
			for line in tempfile:
				line1=line2
				line2=line.rstrip('\n')
				if ">" in line1:
					seqid = line1[1:]
					if revcom == True: seqread = revcom(line2)
					else: seqread = line2
					if seqid not in self.ids:
						aln2 = seqread.find(searchseed)
						self.reads.append(seqread)
						self.ids.append(seqid)
						self.aligns.append(aln1-aln2)
						self.ends.append(aln1-aln2+len(seqread))
						txfound += 1
		return txfound

	def load(self, savefile):
		self.reads = []
		with open(savefile,'r') as readfile:
			intermediary = readfile.readline().rstrip()
			self.parampack = json.loads(intermediary)
			intermediary = readfile.readline().rstrip() 
			self.transcript = json.loads(intermediary)
			intermediary = readfile.readline().rstrip() 
			self.aligns = json.loads(intermediary)
			intermediary = readfile.readline().rstrip() 
			self.ends = json.loads(intermediary)
			intermediary = readfile.readline().rstrip() 
			self.ids = json.loads(intermediary)
			intermediary = readfile.readline().rstrip() 
			self.reads = json.loads(intermediary)
		self.name,self.mindepth,self.kbegin,self.kwidth,self.kend,self.status = self.parampack

	def save(self,savefile):
		self.parampack = [self.name,self.mindepth,self.kbegin,self.kwidth,self.kend,self.status]
		with open(savefile,'w') as writefile:
			json.dump(self.parampack,writefile)
			writefile.write('\n')
			json.dump(self.transcript,writefile)
			writefile.write('\n')
			json.dump(self.aligns,writefile)
			writefile.write('\n')
			json.dump(self.ends,writefile)
			writefile.write('\n')
			json.dump(self.ids,writefile)
			writefile.write('\n')
			json.dump(self.reads,writefile)

	def seed(self,seed,tname):
		self.mindepth = 5
		self.kbegin = 1
		self.kwidth = 20
		self.kend = 21
		self.name = tname
		self.transcript = seed

	def __init__(self):
		self.mindepth = 1
		self.kbegin = 1
		self.kwidth = 20
		self.kend = 21
		self.name = ''
		self.transcript = ''
		self.status = 1
		self.reads = []
		self.ids = []
		self.aligns = []
		self.ends = []

def get_param(prompt_string,screen):
	screen.clear()
	screen.border(0)
	screen.addstr(2, 2, prompt_string)
	screen.refresh()
	input = screen.getstr(10, 10, 60)
	return input

def revcom(inputstr):
	reverse_complement = '' # A reverse compliment is an insult :D
	for n in range(0,len(inputstr)):
		singlechar = inputstr[n:n+1]
		if singlechar == 'A': reverse_complement = 'T' + reverse_complement
		elif singlechar == 'C': reverse_complement = 'G' + reverse_complement
		elif singlechar == 'G': reverse_complement = 'C' + reverse_complement
		elif singlechar == 'T': reverse_complement = 'A' + reverse_complement
		else: reverse_complement = 'N' + reverse_complement
	return reverse_complement

def main(screen,args):
	#								0 WHITE BLACK (n, UNSELECTED)
	curses.init_pair(1,curses.COLOR_GREEN,curses.COLOR_BLACK) #	1 GREEN	BLACK (A, UNSELECTED)
	curses.init_pair(2,curses.COLOR_RED,curses.COLOR_BLACK) #	2 RED	BLACK (C, UNSELECTED)
	curses.init_pair(3,curses.COLOR_MAGENTA,curses.COLOR_BLACK) #	3 MAG	BLACK (G, UNSELECTED)
	curses.init_pair(4,curses.COLOR_CYAN,curses.COLOR_BLACK) #	4 CYAN	BLACK (T, UNSELECTED)
	curses.init_pair(5,curses.COLOR_WHITE,curses.COLOR_BLUE) #	5 WHITE BLUE (n, SELECTED)
	curses.init_pair(6,curses.COLOR_GREEN,curses.COLOR_BLUE) #	6 GREEN BLUE (A, SELECTED)
	curses.init_pair(7,curses.COLOR_RED,curses.COLOR_BLUE) #	7 RED   BLUE (C, SELECTED)
	curses.init_pair(8,curses.COLOR_MAGENTA,curses.COLOR_BLUE) #	8 MAG   BLUE (G, SELECTED)
	curses.init_pair(9,curses.COLOR_CYAN,curses.COLOR_BLUE) #	9 CYAN  BLUE (T, SELECTED)
	curses.init_pair(10,curses.COLOR_WHITE,curses.COLOR_YELLOW) #	10 WHITE YELLOW (n, HIGHLIGHTED)
	curses.init_pair(11,curses.COLOR_GREEN,curses.COLOR_YELLOW) #	11 GREEN YELLOW (A, HIGHLIGHTED)
	curses.init_pair(12,curses.COLOR_RED,curses.COLOR_YELLOW) #	12 RED   YELLOW (C, HIGHLIGHTED)
	curses.init_pair(13,curses.COLOR_MAGENTA,curses.COLOR_YELLOW) #	13 MAG   YELLOW (G, HIGHLIGHTED)
	curses.init_pair(14,curses.COLOR_CYAN,curses.COLOR_YELLOW) #	14 CYAN  YELLOW (T, HIGHLIGHTED)
	# window init complete
	strandlist = []
	extendonly = -1
	#File init
	inputlist = [args.onepath,args.twopath,args.threepath,args.fourpath]
	while 1:
		screen.clear()
		screen.border(0)
		screen.addstr(2, 2, "Please enter a number...")
		screen.addstr(4, 4, "0 - Set file inputs")
		screen.addstr(5, 4, "1 - Initialize Transcript From Seed")
		screen.addstr(6, 4, "2 - Extend Transcript")
		screen.addstr(7, 4, "3 - View Transcript")
		screen.addstr(8, 4, "4 - Save File")
		screen.addstr(9, 4, "5 - Load File")
		screen.addstr(10, 4, "6 - Print Data")
		screen.addstr(11, 4, "7 - Set Extend-only")
		screen.addstr(12, 4, "q - Exit")
		screen.refresh()

		x = screen.getch()
		if x == ord('q'):
			break
		elif x == ord('0'): #file options
			while 1:
				screen.clear()
				screen.border(0)
				screen.addstr(2, 2, "Please enter a number...")
				screen.addstr(4, 4, "1 - " + inputlist[0])
				screen.addstr(5, 4, "2 - " + inputlist[1])
				screen.addstr(6, 4, "3 - " + inputlist[2])
				screen.addstr(7, 4, "4 - " + inputlist[3])
				screen.addstr(8, 4, "q - Exit")
				screen.refresh()
				z = screen.getch()
				if z== ord('q'):
					break
				elif z==ord('1'):
					if inputlist[0] == args.onepath: inputlist[0] = args.nullpath
					else: inputlist[0] = args.onepath
				elif z==ord('2'):
					if inputlist[1] == args.twopath: inputlist[1] = args.nullpath
					else: inputlist[1] = args.twopath
				elif z==ord('3'):
					if inputlist[2] == args.threepath: inputlist[2] = args.nullpath
					else: inputlist[2] = args.threepath
				elif z==ord('4'):
					if inputlist[3] == args.fourpath: inputlist[3] = args.nullpath
					else: inputlist[3] = args.fourpath
		elif x == ord('1'): #seed
			screen.clear()
			screen.border(0)
			screen.addstr(2, 2, "Initializing...")
			screen.refresh()
			strandlist = []
			txnum = 0
			line1 = ''
			line2 = ''
			with open(args.seedpath,'r') as seedfile:
				for line in seedfile:
					line1=line2
					line2=line.rstrip('\n')
					if ">" in line1:
						strandlist.append(BLASToys())
						txnum += 1
						strandlist[-1].seed(line2,str(txnum)+'_'+line1[1:])
		elif x == ord('2'):
#			runlist = get_param("Run the following:" ,screen)
			screen.clear()
			screen.border(0)
			screen.addstr(2, 2, "Extending...")
			screen.refresh()
			runset = []
			if extendonly == -1: runset = range(0,len(strandlist))
			else: runset.append(extendonly)
			for filepath in inputlist:
				for p in runset:
					txfound = 0
					screen.addstr(3+p%20, 2, str(p)+'a')
					screen.refresh()
					txfound += strandlist[p].transextend(filepath,args.temppath,False)
					screen.addstr(3+p%20, 2, str(p)+'b')
					screen.refresh()
					txfound += strandlist[p].transextend(filepath,args.temppath,True)
					screen.addstr(3+p%20,10, 'Tx found: ' + str(txfound))
					screen.refresh()
					if txfound > 0:
						if min(strandlist[p].aligns) < 0 :
							padsize = abs(min(strandlist[p].aligns))
							strandlist[p].transcript = 'n'*padsize+strandlist[p].transcript
							strandlist[p].kbegin += padsize
							strandlist[p].kend += padsize
							for n in range(0,len(strandlist[p].reads)):
								strandlist[p].aligns[n] += padsize
								strandlist[p].ends[n] += padsize
						if max(strandlist[p].ends) > len(strandlist[p].transcript):
							padsize = max(strandlist[p].ends) - len(strandlist[p].transcript)
							strandlist[p].transcript = strandlist[p].transcript)+'n'*padsize
			screen.getch()
		elif x == ord('3'):
			screen.clear()
			screen.border(0)
			screen.addstr(2, 2, "Rendering...")
			screen.refresh()
			stranddex = 0
			while stranddex != -1:
				stranddex = transrender(screen,strandlist,stranddex)
		elif x == ord('4'): #save
			for p in range(0,len(strandlist)):
				outfilename = args.savepath + '/' + str(p+1000) + '.json'
				strandlist[p].save(outfilename)
		elif x == ord('5'): #load
			#list files
			loadlist = []
			strandlist = []
			(dirpath, _, filenames) = os.walk(args.loadpath).next()
			for filename in filenames:
				if ".json" in filename: loadlist.append(int(filename[:-5]))
			loadlist.sort()
			#iterate over files found
			for p in range(0,len(loadlist)):
				infilename = args.loadpath+'/'+str(loadlist[p])+'.json'
				strandlist.append(BLASToys())
				strandlist[-1].load(infilename)
		elif x == ord('6'): #print
			with open(args.printpath,'w') as printfile:
				for p in range(0,len(strandlist)):
					printfile.write('>'+ strandlist[p].name + '\n')
					printfile.write(strandlist[p].transcript + '\n')
		elif x == ord('7'): #set extendonly
			extendonly = 0			
	curses.endwin()
	exit()

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('-1', '--one', dest='onepath')
	parser.add_argument('-2', '--two', dest='twopath')
	parser.add_argument('-3', '--three', dest='threepath')
	parser.add_argument('-4', '--four', dest='fourpath')
	parser.add_argument('-n', '--null', dest='nullpath')
	parser.add_argument('-x', '--seed', dest='seedpath')
	parser.add_argument('-s', '--save', dest='savepath')
	parser.add_argument('-l', '--load', dest='loadpath')
	parser.add_argument('-p', '--print', dest='printpath')
	parser.add_argument('-t', '--temp', dest='temppath')
	args = parser.parse_args()
	curses.wrapper(main,args)
