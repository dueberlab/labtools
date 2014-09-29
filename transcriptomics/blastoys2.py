#!/usr/bin/env python
#Written by ZNR; 2014/09/07
import argparse
import curses
import copy
import os
import json
import subprocess
import sys

#Shared stuff here
class BLASToys:
	def __init__(self):
		self.mindepth = 1
		self.minsplit = 1
		self.kbegin = 1
		self.kwidth = 31
		self.kend = 32
		self.name = ''
		self.transcript = ''
		self.status = 1
		self.list_read_sequences = []
		self.list_read_ids = []
		self.list_read_aligns = []
		self.list_read_ends = []
		self.coverage = []

	def deleteReads(self,list_indices):
		self.list_read_sequences = deleteByIndices(self.list_read_sequences, list_indices)
		self.list_read_ids = deleteByIndices(self.list_read_ids, list_indices)
		self.list_read_aligns = deleteByIndices(self.list_read_aligns, list_indices)
		self.list_read_ends = deleteByIndices(self.list_read_ends, list_indices)

	def findReads(self, inputlist, iovar_temppath, enablerc):
		searchseed = self.transcript[self.kbegin:self.kend]
		aln1 = self.transcript.find(searchseed)
		local_padsize = 0
		local_txfound = 0
		rcbool = [False]
		if enablerc == True: rcbool.append(True)
		# must have single line, non-wrapping fastas to work - no n60 splitting; otherwise, can't use grep for max speed
		for local_filepath in inputlist:
			for rcset in rcbool:
				with open(iovar_temppath + '/temp','w') as tempfile:
					if rcset == True:
						grep = subprocess.Popen(['grep', '-F', '-B1', reverseComplement(searchseed), filepath],stdout=tempfile)
					else:
						grep = subprocess.Popen(['grep', '-F', '-B1', searchseed, filepath],stdout=tempfile)
					grep.wait()
					tempfile.flush()
				with open(iovar_temppath + '/temp','r') as tempfile:
					line1=''
					line2=''
					for line in tempfile:
						line1=line2
						line2=line.rstrip('\n')
						if ">" in line1:
							seqid = line1[1:]
							if rcset == True: seqread = reverseComplement(line2)
							else: seqread = line2
							if seqid not in self.list_read_ids:
								aln2 = seqread.find(searchseed)
								self.list_read_sequences.append(seqread.upper())
								self.list_read_ids.append(seqid)
								self.list_read_aligns.append(aln1-aln2)
								self.list_read_ends.append(aln1-aln2+len(seqread))
								local_txfound += 1
		if local_txfound > 0:
			self.coverage[self.kbegin:self.kend] = [local_txfound] * self.kwidth
			if min(self.list_read_aligns) < 0 :
				local_padsize = abs(min(self.list_read_aligns))
				self.transcript = 'n'*local_padsize+self.transcript
				self.coverage = [0]*local_padsize+self.coverage
				self.kbegin += local_padsize
				self.kend += local_padsize
				for n in range(0,len(self.list_read_sequences)):
					self.list_read_aligns[n] += local_padsize
					self.list_read_ends[n] += local_padsize
			if max(self.list_read_ends) > len(self.transcript):
				local_padsize = max(self.list_read_ends) - len(self.transcript)
				self.transcript = self.transcript+'n'*local_padsize
				self.coverage = self.coverage + [0] * local_padsize
		else: self.coverage[self.kbegin:self.kend] = [-1] * self.kwidth
		return local_txfound
		
	def generateReference(self):
		newtranscript = ''
		newcoverage = []
		newrefs = []
		for n in range(0,len(self.transcript)):
			ctA=0
			ctC=0
			ctG=0
			ctT=0
			readsA = []
			readsC = []
			readsG = []
			readsT = []
			readsN = []
			for o in range(0,len(self.list_read_sequences)):
				yread = self.list_read_sequences[o]
				yalign = self.list_read_aligns[o]
				yend = self.list_read_ends[o]
				if yalign<=n<yend:
					basechar = yread[n-yalign:n-yalign+1]
					if basechar in 'A':
						ctA +=1
						readsA.append(o)
					elif basechar in 'C':
						ctC +=1
						readsC.append(o)
					elif basechar in 'G':
						ctG +=1
						readsG.append(o)
					elif basechar in 'T':
						ctT +=1
						readsT.append(o)
					else:
						readsN.append(o)
			maxN = max([ctA,ctC,ctG,ctT])
			if maxN >= self.mindepth:
				if maxN == ctA:
					newtranscript += 'A'
					if max([ctC,ctG,ctT]) >= self.minsplit:
						newrefs.append(copy.deepcopy(self))
						newrefs[-1].name = self.name+"_notA"
						newrefs[-1].deleteReads(readsA)
						self.deleteReads([readsC,readsG,readsT,readsN])
				elif maxN == ctC:
					newtranscript += 'C'
					if max([ctA,ctG,ctT]) >= self.minsplit:
						newrefs.append(copy.deepcopy(self))
						newrefs[-1].name = self.name+"_notC"
						newrefs[-1].deleteReads(readsC)
						self.deleteReads([readsA,readsG,readsT,readsN])
				elif maxN == ctG:
					newtranscript += 'G'
					if max([ctA,ctC,ctT]) >= self.minsplit:
						newrefs.append(copy.deepcopy(self))
						newrefs[-1].name = self.name+"_notG"
						newrefs[-1].deleteReads(readsG)
						self.deleteReads([readsA,readsC,readsT,readsN])
				elif maxN == ctT:
					newtranscript += 'T'
					if max([ctA,ctC,ctG]) >= self.minsplit:
						newrefs.append(copy.deepcopy(self))
						newrefs[-1].name = self.name+"_notT"
						newrefs[-1].deleteReads(readsT)
						self.deleteReads([readsA,readsC,readsG,readsN])
			else:
				newtranscript += 'n'
		self.transcript=newtranscript
		return newrefs

	def kmerScan(self):
		list_available = []
		for loop_index in range(0,len(self.coverage)):
			if [0]*self.kwidth == self.coverage[loop_index:loop_index+self.kwidth] and 'n' not in self.transcript[loop_index:loop_index+self.kwidth]: list_available.append(loop_index)
		if list_available[0:1]:
			self.kbegin = list_available[0:1]
			self.kend = self.kbegin+self.kwidth
			return 1
		else
			return 0
		
		
	def loadJSON(self, savefile):
		self.list_read_sequences = []
		with open(savefile,'r') as readfile:
			self.parampack = json.loads(readfile.readline().rstrip())
			self.transcript = json.loads(readfile.readline().rstrip())
			self.list_read_aligns = json.loads(readfile.readline().rstrip())
			self.list_read_ends = json.loads(readfile.readline().rstrip())
			self.list_read_ids = json.loads(readfile.readline().rstrip())
			self.coverage = json.loads(readfile.readline().rstrip())
			self.list_read_sequences = json.loads(readfile.readline().rstrip())
		self.name,self.mindepth,self.minsplit,self.kbegin,self.kwidth,self.kend,self.status = self.parampack

	def saveJSON(self,savefile):
		self.parampack = [self.name,self.mindepth,self.minsplit,self.kbegin,self.kwidth,self.kend,self.status]
		with open(savefile,'w') as writefile:
			json.dump(self.parampack,writefile)
			writefile.write('\n')
			json.dump(self.transcript,writefile)
			writefile.write('\n')
			json.dump(self.list_read_aligns,writefile)
			writefile.write('\n')
			json.dump(self.list_read_ends,writefile)
			writefile.write('\n')
			json.dump(self.list_read_ids,writefile)
			writefile.write('\n')
			json.dump(self.coverage,writefile)
			writefile.write('\n')
			json.dump(self.list_read_sequences,writefile)

	def seed(self,seed,tname):
		self.mindepth = 3
		self.minsplit = 3
		self.kwidth = 31
		self.name = tname
		self.transcript = seed
		self.kbegin = int(math.floor(len(self.transcript)/2))
		self.kend = self.kbegin+self.kwidth
		self.coverage = [0 for range(0,len(self.transcript))]
	
def deleteByIndices(lst, indices):
	# from http://stackoverflow.com/questions/497426/deleting-multiple-elements-from-a-list
	return [ lst[i] for i in range(0,len(lst)) if i not in set(indices) ]

def reverseComplement(inputstr):
	reverse_complement = '' # A reverse compliment is an insult :D
	for n in range(0,len(inputstr)):
		singlechar = inputstr[n:n+1]
		if singlechar == 'A': reverse_complement = 'T' + reverse_complement
		elif singlechar == 'C': reverse_complement = 'G' + reverse_complement
		elif singlechar == 'G': reverse_complement = 'C' + reverse_complement
		elif singlechar == 'T': reverse_complement = 'A' + reverse_complement
		else: reverse_complement = 'N' + reverse_complement
	return reverse_complement

def loadTranscripts(loadpath):
	#list files
	loadlist = []
	list_all_strands = []
	(dirpath, _, filenames) = os.walk(loadpath).next()
	for filename in filenames:
		if ".json" in filename: loadlist.append(int(filename[:-5]))
		loadlist.sort()
		#iterate over files found
		for loop_p in range(0,len(loadlist)):
			infilename = loadpath+'/'+str(loadlist[loop_p])+'.json'
			list_all_strands.append(BLASToys())
			list_all_strands[-1].loadJSON(infilename)
	return list_all_strands

def saveTranscripts(list_all_strands,savepath):
	for loop_p in range(0,len(list_all_strands)):
		outfilename = savepath + '/' + str(loop_p+1000000) + '.json'
		list_all_strands[p].saveJSON(outfilename)

def printTranscripts(list_all_strands,printpath):
	with open(args.printpath,'w') as printfile:
		for loop_p in range(0,len(list_all_strands)):
			printfile.write('>'+ list_all_strands[loop_p].name + '\n')
			printfile.write(list_all_strands[loop_p].transcript + '\n')

def seedTranscripts(seedpath):
	list_all_strands = []
	txnum = 0
	line1 = ''
	line2 = ''
	with open(seedpath,'r') as seedfile:
		for line in seedfile:
			line1=line2
			line2=line.rstrip('\n')
			if ">" in line1:
				list_all_strands.append(BLASToys())
				txnum += 1
				list_all_strands[-1].seed(line2,str(txnum)+'_'+line1[1:])
	return list_all_strands

def getParam(prompt_string,screen):
	screen.clear()
	screen.border(0)
	screen.addstr(2, 2, prompt_string)
	screen.refresh()
	input = screen.getstr(10, 10, 60)
	return input

def renderTranscript(screen,list_all_strands,index_strandlist):
	#Reset cursor and view position
	subindex_base = 0
	subindex_reads = 0
	render_xview = 0
	render_yview = 0
	#Fix graphics boundaries
	leftbound = 1
	upperbound = 6
	lowerbound=curses.LINES-1
	rightbound=curses.COLS-1
	#Misc vars
	element_strand = list_all_strands[index_strandlist]
	sellist = []
	while 1:
		#Calculate cursor locations
		render_xcursor = subindex_base-render_xview
		render_ycursor = subindex_reads-render_yview
		numreads = len(element_strand.list_read_sequences)
		element_strand.kend = element_strand.kbegin + element_strand.kwidth
		# Refresh screen and labels
		screen.clear()
		screen.border(0)
		screen.addstr(1, 1, "Strand ID:"+str(index_strandlist)+"of"+str(len(list_all_strands)-1))
		screen.addstr(1, 20, "Min depth:"+str(element_strand.mindepth))
		screen.addstr(1, 35, "kmer width:"+str(element_strand.kwidth))
		screen.addstr(2, 1, "q to quit, wasd to move view (shift to go faster), op to switch strand, arrows to set kmer bounds, -=[] to move selector (shift to go faster)")
		screen.addstr(3, 1, "r to generate reference,  m to switch min depth, n to copy to new transcript, x to remove read, shift-Z to delete transcript")
		# Blit
		for xblit in range(leftbound,rightbound):
			xindex = xblit-leftbound+render_xview
			#individual strand status

			#transcript
			if xindex == subindex_base: rcolor = 10
			elif xindex in range(element_strand.kbegin,element_strand.kend): rcolor = 5
			else: rcolor = 0
			basechar = element_strand.transcript[xindex:xindex+1]
			yblit = 5
			if 'A' in basechar: screen.addstr(yblit,xblit,'A',curses.color_pair(rcolor+1))
			elif 'C' in basechar: screen.addstr(yblit,xblit,'C',curses.color_pair(rcolor+2))
			elif 'G' in basechar: screen.addstr(yblit,xblit,'G',curses.color_pair(rcolor+3))
			elif 'T' in basechar: screen.addstr(yblit,xblit,'T',curses.color_pair(rcolor+4))
			else: screen.addstr(yblit,xblit,'n',curses.color_pair(rcolor))
			#reads
			for yblit in range(upperbound,lowerbound):
				yindex = yblit-upperbound+render_yview
				if yindex < numreads:
					yread = element_strand.list_read_sequences[yindex]
					yalign = element_strand.list_read_aligns[yindex]
					yend = element_strand.list_read_ends[yindex]
					#Row Highlighting
					if yindex == subindex_reads or xindex == subindex_base: rcolor = 10
					elif yindex in sellist: rcolor = 5
					else: rcolor = 0
					#Blit bits
					if xindex<yalign: screen.addstr(yblit,xblit,'.',curses.color_pair(rcolor))
					elif xindex>=yend: screen.addstr(yblit,xblit,'.',curses.color_pair(rcolor))
					elif yalign<=xindex<yend:
						basechar = yread[xindex-yalign:xindex-yalign+1]
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
			index_strandlist = -1
			break
		elif y == ord('o') and index_strandlist != 0:
			index_strandlist += -1
			break
		elif y == ord('p') and index_strandlist+1 < len(list_all_strands):
			index_strandlist += 1
			break
		elif y == ord('-') and subindex_base > 0: subindex_base +=-1	# Navigation options
		elif y == ord('_') and subindex_base > 5: subindex_base +=-5
		elif y == ord('=') and subindex_base < len(element_strand.transcript): subindex_base +=1
		elif y == ord('+') and subindex_base < len(element_strand.transcript)-5: subindex_base +=5
		elif y == ord('[') and subindex_reads > 0: subindex_reads +=-1
		elif y == ord('{') and subindex_reads > 5: subindex_reads +=-5
		elif y == ord(']') and subindex_reads+1 < numreads: subindex_reads +=1
		elif y == ord('}') and subindex_reads+1 < numreads-5: subindex_reads +=5
		elif y == ord('w') and render_yview>0: render_yview +=-1
		elif y == ord('W') and render_yview>5: render_yview +=-5
		elif y == ord('a') and render_xview>0: render_xview +=-1
		elif y == ord('A') and render_xview>5: render_xview +=-5
		elif y == ord('s') and render_yview<(numreads): render_yview +=1
		elif y == ord('S') and render_yview<(numreads)-5: render_yview +=5
		elif y == ord('d') and render_xview+(rightbound-leftbound)<len(element_strand.transcript): render_xview +=1
		elif y == ord('D') and render_xview+(rightbound-leftbound)<len(element_strand.transcript)-5: render_xview +=5
		elif y == curses.KEY_LEFT and element_strand.kbegin>0: element_strand.kbegin +=-1
		elif y == curses.KEY_RIGHT and element_strand.kwidth+element_strand.kbegin<len(element_strand.transcript): element_strand.kbegin  +=1
		elif y == curses.KEY_UP and element_strand.kwidth>0: element_strand.kwidth +=-1
		elif y == curses.KEY_DOWN and element_strand.kwidth+element_strand.kbegin<len(element_strand.transcript): element_strand.kwidth +=1
		elif y == ord(' '):	# Misc Options
			if subindex_reads in sellist: sellist.remove(subindex_reads)
			else: sellist.append(subindex_reads)
		elif y == ord('n'):
			list_all_strands.append(copy.deepcopy(element_strand))
			list_all_strands[-1].name = str(len(list_all_strands)-1)+"_"+element_strand.name
		elif y == ord('Z'):
			list_all_strands.pop(index_strandlist)
			index_strandlist += -1
			break
		elif y == ord('x'):
			element_strand.deleteReads(subindex_reads)
			sellist.pop(subindex_reads)
			subindex_reads +=-1
			if subindex_reads < 0: subindex_reads=0
		elif y == ord('r'):
			list_all_strands.append(element_strand.generateReference())
		elif y == ord('m'):
			depths = [1,2,3,4,5,10]
			element_strand.mindepth = depths[depths.index(element_strand.mindepth)-1]
	return index_strandlist

def automatedMain(args):
	loop_runs = 0
	inputlist = args.iovar_inputpaths.split(',')
	if hasattr(args, 'loadpath'): list_all_strands = loadTranscripts(args.loadpath)
	elif hasattr(args, 'seedpath'): list_all_strands = seedTranscripts(args.seedpath)
	if hasattr(args, 'maxiter'):
		while loop_runs < args.maxiter:
			for element_strand in list_all_strands:
				if element_strand.status > 0 and element_strand.kmerScan() == 0:  element_strand.status = 0
				if element_strand.status > 0 and element_strand.findReads(inputlist,args.iovar_temppath,True) == 0: element_strand.status += 1
				if element_strand.status > 10: element_strand.status = 0
				if element_strand.status > 0: list_all_strands.append(element_strand.generateReference())
			loop_runs += 1
	if hasattr(args, 'savepath'): saveTranscripts(list_all_strands,args.savepath):
	if hasattr(args, 'printpath'): printTranscripts(list_all_strands,args.printpath):
	exit()


def interactiveMain(screen,args):
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
	list_all_strands = []
	extendonly = -1
	#File init
	inputlist = args.iovar_inputpaths.split(',')
	while 1:
		screen.clear()
		screen.border(0)
		screen.addstr(2, 2, "Please enter a number...")
		screen.addstr(4, 4, "0 - Set file inputs")
		screen.addstr(5, 4, "1 - Seed or Load")
		screen.addstr(6, 4, "2 - Extend Transcript")
		screen.addstr(7, 4, "3 - View Transcript")
		screen.addstr(8, 4, "4 - Save File")
		screen.addstr(9, 4, "5 - Print Data")
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
				for printlist in range(0,len(inputlist)):
					screen.addstr(4, 4, str(printlist)+" - " + inputlist[printlist])
				screen.addstr(30, 4, "q - Exit")
				screen.refresh()
				z = screen.getch()
				if z== ord('q'):
					break

		elif x == ord('1'): #seed or load
			screen.clear()
			screen.border(0)
			screen.addstr(2, 2, "Initializing...")
			screen.refresh()
			if hasattr(args, 'loadpath'): list_all_strands = loadTranscripts(args.loadpath)
			elif hasattr(args, 'seedpath'): list_all_strands = seedTranscripts(args.seedpath)
			
		elif x == ord('2'):
#			runlist = getParam("Run the following:" ,screen)
			screen.clear()
			screen.border(0)
			screen.addstr(2, 2, "Extending...")
			screen.refresh()
			for element_strand in list_all_strands:
				if element_strand.status > 0: local_txfound = element_strand.findReads(inputlist,args.iovar_temppath,True)
			
		elif x == ord('3'):
			screen.clear()
			screen.border(0)
			screen.addstr(2, 2, "Rendering...")
			screen.refresh()
			index_strandlist = 0
			while index_strandlist != -1:
				index_strandlist = renderTranscript(screen,list_all_strands,index_strandlist)
		
		elif x == ord('4'):
			screen.clear()
			screen.border(0)
			screen.addstr(2, 2, "Saving...")
			screen.refresh()
			if hasattr(args, 'savepath'): saveTranscripts(list_all_strands,savepath):
		
		elif x == ord('5'):
			screen.clear()
			screen.border(0)
			screen.addstr(2, 2, "Printing...")
			screen.refresh()
			if hasattr(args, 'printpath'): printTranscripts(list_all_strands,printpath):
		
		elif x == ord('7'): #set extendonly
			extendonly = 0			

	curses.endwin()
	exit()

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('-a', '--auto', action='store_true')
	parser.add_argument('-q', '--quiet', action='store_true')
	parser.add_argument('-i', '--inputlist', dest='inputlist')
	parser.add_argument('-m', '--maxiter', dest='maxiter')
	parser.add_argument('-x', '--seed', dest='seedpath')
	parser.add_argument('-s', '--save', dest='savepath')
	parser.add_argument('-i', '--inputpaths', dest='iovar_inputpaths')
	parser.add_argument('-l', '--load', dest='loadpath')
	parser.add_argument('-p', '--print', dest='printpath')
	parser.add_argument('-t', '--temp', dest='iovar_temppath')
	args = parser.parse_args()
	if args.auto:
		automatedMain(args)
	else:
		curses.wrapper(interactiveMain,args)
