import struct

def stringPad(string, character, length, direction="left"):
	string = str(string)
	if direction == "right":
		return string.ljust(length, character)
	else:
		return string.rjust(length, character)

def to8BitArray(val):
	return bytearray(struct.pack('B', val))

def to16BitLEArray(val):
	return bytearray(struct.pack('<H', val))

def to32BitLEArray(val):
	return bytearray(struct.pack('<I', val))

def to16BitBEArray(val):
	return bytearray(struct.pack('>H', val))

def to32BitBEArray(val):
	return bytearray(struct.pack('>I', val))
