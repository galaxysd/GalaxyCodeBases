import random
import re
import struct

maxUInt = 4294967295

def randomUInt32(*args):
	if len(args) == 0:
		return random.randint(0, maxUInt)
	elif len(args) == 1:
		maximum = args[0]
		minimum = 0
	else:
		minimum = args[0]
		maximum = args[1]
	resultRange = (maximum + 1) - minimum
	factor = resultRange / maxUInt
	return int((random.randint(0, maxUInt) * factor) + minimum)

def stringPad(string, character, length, direction="left"):
	string = str(string)
	pad = character.join(arrayFill([], max(1, length - len(string) + 1), "")) 
	if direction == "right":
		return (string + pad)[:length]
	else:
		return (pad + string)[-1 * length:]

def arrayFill(arr, length, value):
	while length > 0:
		arr.append(value)
		length -= 1
	return arr

def byteArrayToString(arr):
	string = ""
	for i in range(0,len(arr)):
		string += chr(arr[i])
	return string

def bufferToByteArray(buf, start=0, length=None):
	if not length:
		length = len(buf)
	end = start + length
	arr = bytearray()
	for i in range(start,end):
		arr.append(buf[i])
	return arr

def to8BitArray(val):
	arr = bytearray()
	bytes = struct.pack('B',val)
	for i in range(0,len(bytes)):
		arr.append(bytes[i])
	return arr

def to16BitLEArray(val):
	arr = bytearray()
	bytes = struct.pack('<H',val)
	for i in range(0,len(bytes)):
		arr.append(bytes[i])
	return arr

def to32BitLEArray(val):
	arr = bytearray()
	bytes = struct.pack('<I',val)
	for i in range(0,len(bytes)):
		arr.append(bytes[i])
	return arr

def to16BitBEArray(val):
	arr = bytearray()
	bytes = struct.pack('>H',val)
	for i in range(0,len(bytes)):
		arr.append(bytes[i])
	return arr

def to32BitBEArray(val):
	arr = bytearray()
	bytes = struct.pack('>I',val)
	for i in range(0,len(bytes)):
		arr.append(bytes[i])
	return arr
