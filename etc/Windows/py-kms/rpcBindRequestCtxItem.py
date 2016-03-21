import functions
import struct
import uuid

class rpcBindRequestCtxItem:
	dataLength = 44

	def __init__(self, data):
		self.contextId = struct.unpack_from('<H', str(data), 0)[0];
		self.numTransItems = struct.unpack_from('<H', str(data), 2)[0];
		self.abstractSyntaxUUID = uuid.UUID(bytes_le=str(data[4:4 + 16]))
		self.abstractSyntaxVersionMajor = struct.unpack_from('<H', str(data), 20)[0]
		self.abstractSyntaxVersionMinor = struct.unpack_from('<H', str(data), 22)[0]
		self.transferSyntaxUUID = uuid.UUID(bytes_le=str(data[24:24 + 16]))
		self.transferSyntaxVersion = struct.unpack_from('<I', str(data), 40)[0]

	def toDictionary(self):
		dictionary = {}
		dictionary['contextId'] = self.contextId
		dictionary['numTransItems'] = self.numTransItems
		dictionary['abstractSyntaxUUID'] = self.abstractSyntaxUUID
		dictionary['abstractSyntaxVersionMajor'] = self.abstractSyntaxVersionMajor
		dictionary['abstractSyntaxVersionMinor'] = self.abstractSyntaxVersionMinor
		dictionary['transferSyntaxUUID'] = self.transferSyntaxUUID
		dictionary['transferSyntaxVersion'] = self.transferSyntaxVersion
		return dictionary

	def toByteArray(self):
		bytes = bytearray()
		bytes.extend(functions.to16BitLEArray(self.contextId))
		bytes.extend(functions.to16BitLEArray(self.numTransItems))
		bytes.extend(self.abstractSyntaxUUID.bytes_le)
		bytes.extend(functions.to16BitLEArray(self.abstractSyntaxVersionMajor))
		bytes.extend(functions.to16BitLEArray(self.abstractSyntaxVersionMinor))
		bytes.extend(self.transferSyntaxUUID.bytes_le)
		bytes.extend(functions.to32BitLEArray(self.transferSyntaxVersion))
		return bytes
