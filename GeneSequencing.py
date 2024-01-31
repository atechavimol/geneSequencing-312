#!/usr/bin/python3

from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
	from PyQt4.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import math
import time
import random

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1

class GeneSequencing:

	def __init__( self ):
		pass
	
	'''
	Time Complexity: O(nm) or O(kn)
	Space Complexity: O(nm) or O(kn)
	'''
	def align( self, seq1, seq2, banded, align_length):
		self.banded = banded
		self.MaxCharactersToAlign = align_length

		# depending on the whether banded or not, call the appropriate alignment function
		if banded:
			sol, seq1, seq2 = self.find_alignment_banded(seq1[:align_length], seq2[:align_length])
		else:
			sol, seq1, seq2 = self.find_alignment_unbanded(seq1[:align_length], seq2[:align_length])

		# returns best cost and sequence alignment strings
		return {'align_cost': sol, 'seqi_first100':seq1, 'seqj_first100':seq2}

	'''
		Time Complexity: O(nm)
		Space Complexity: O(nm)
	'''
	def find_alignment_unbanded(self, seq1, seq2):
		self.dict = {}

		# check if they are the same seq
		if seq1 == seq2:
			return len(seq1) * MATCH, seq1, seq2

		seq1 = "*" + seq1
		seq2 = "*" + seq2

		for i in range(0, len(seq1)):
			for j in range(0, len(seq2)):
				# initializes top row
				if i == 0:
					if j == 0:
						self.dict[(0, j)] = Cell(INDEL * j, None)
					else:
						self.dict[(0, j)] = Cell(INDEL * j, (0, j -1))
				# initializes leftmost column
				elif j == 0:
					if i == 0:
						self.dict[(i, 0)] = Cell(INDEL * i, None)
					else:
						self.dict[(i, 0)] = Cell(INDEL * i, (0, i - 1))
				else:
					# calculates the cost for if seq1[i] == seq2[j] matches or not
					diff = SUB
					if seq1[i] == seq2[j]:
						diff = MATCH

					# finds whether the Cell top, diagonal or left has the smallest cost
					min_coord, min_val = self.get_min(i, j, diff)

					# creates cell with lowest cost pointing back to where it came from
					self.dict[(i,j)] = Cell(min_val, min_coord)

		#gets optimal alignment based on the Cell's prev coordinate
		seq1_alignment, seq2_alignment = self.traceback(seq1, seq2)

		# returns the total cost and alignment strings
		return self.dict[(len(seq1) - 1, len(seq2) - 1)].val, seq1_alignment[:100], seq2_alignment[:100]

	'''
		Time Complexity: O(kn)
		Space Complexity: O(kn)
	'''
	def find_alignment_banded(self, seq1, seq2):
		self.dict = {}
		seq1 = "*" + seq1
		seq2 = "*" + seq2
		d=3
		k=7

		# returns infinity and 'No Alignment Possible' if the size difference is greater than the band's scope
		if abs(len(seq1) - len(seq2)) > d:
			return math.inf, "No Alignment Possible", "No Alignment Possible"

		# initializes the top row and leftmost column
		for i in range(d + 1):
			if i == 0:
				self.dict[(i, 0)] = Cell(i * INDEL, None)
				self.dict[(0, i)] = Cell(i * INDEL, None)
			else:
				self.dict[(i, 0)] = Cell(i * INDEL, (i-1, 0))
				self.dict[(0, i)] = Cell(i * INDEL, (0, i-1))


		for i in range(1, len(seq1)):
			for x in range(i-d, i+d + 1):
				# checks if x is a valid index in the table
				if x > 0 and x < len(seq2):
					# calculates the cost for if seq1[i] == seq2[j] matches or not
					diff = SUB
					if seq1[i] == seq2[x]:
						diff = MATCH

					# finds whether the Cell top, diagonal or left has the smallest cost
					min_coord, min_val = self.get_min(i, x, diff)
					# creates cell with lowest cost pointing back to where it came from
					self.dict[(i,x)] = Cell(min_val, min_coord)

		#gets optimal alignment based on the Cell's prev coordinate
		seq1_alignment, seq2_alignment = self.traceback(seq1, seq2)

		# returns the total cost and alignment strings


		return self.dict[(len(seq1) - 1, len(seq2) - 1)].val, seq1_alignment[:100], seq2_alignment[:100]

	'''
		Time Complexity: O(1)
		Space Complexity: O(1)
	'''
	def get_min(self, i, j, diff):
		# initializes left and top to infinity (assuming that the cell(i,j) doesn't have a cell to the left or above it,
		# therefore having a cost of infinity meaning unreachable
		left_val = math.inf
		top_val = math.inf

		# if there exists a cell to the left of (i,j), get it's cost val
		if self.dict.get((i-1, j)) is not None:
			left_val = self.dict[(i-1, j)].val + INDEL

		# if there exists a cell on top of (i,j), get it's cost val
		if self.dict.get((i, j-1)) is not None:
			top_val = self.dict[(i, j-1)].val + INDEL

		# get the cost of the diagonal cell
		diag_val = self.dict[(i-1, j-1)].val + diff

		# with a tie breaking precedence of left, top, and diag, return the smallest cost and coordinate of that smallest cost
		if left_val <= top_val and left_val <= diag_val:
			return (i-1, j), left_val
		elif top_val <= diag_val:
			return (i, j-1), top_val
		else:
			return (i-1, j-1), diag_val

	'''
		Time Complexity: O(m + n)
		Space Complexity: O(m + n)
	'''
	def traceback(self, seq1, seq2):
		seq1_alignment = ""
		seq2_alignment = ""

		# position at bottom right-hand corner
		curr_pos = (len(seq1) - 1, len(seq2) - 1)

		# position of curr_pos's previous coordinate
		prev_pos = self.dict[(len(seq1) - 1, len(seq2) - 1)].prev

		while prev_pos is not None:

			# if prev is diagonal to curr, add letters from seq1 and seq1 into alignment strings
			if prev_pos == (curr_pos[0] - 1, curr_pos[1] - 1):
				seq1_alignment = seq1[curr_pos[0]] + seq1_alignment
				seq2_alignment = seq2[curr_pos[1]] + seq2_alignment

			# if prev is on top of curr, insert into seq1 aligment string and add letter from seq2 into alignment string
			elif prev_pos == (curr_pos[0], curr_pos[1] - 1):
				seq1_alignment = "-" + seq1_alignment
				seq2_alignment = seq2[curr_pos[1]] + seq2_alignment

			# if prev is to the left curr, insert into seq2 aligment string and add letter from seq1 into alignment string
			else:
				seq1_alignment = seq1[curr_pos[0]] + seq1_alignment
				seq2_alignment = "-" + seq2_alignment

			# update positions
			curr_pos = prev_pos
			prev_pos = self.dict[curr_pos].prev

		return seq1_alignment, seq2_alignment


class Cell:
	def __init__(self, val, prev):
		self.val = val
		self.prev = prev