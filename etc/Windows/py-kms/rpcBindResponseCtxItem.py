import functions

class rpcBindResponseCtxItem:
	def __init__(self, data):
		self.ackResult = data['ackResult']
		self.ackReason = data['ackReason']
		self.transferSyntax = data['transferSyntax']
		self.syntaxVersion = data['syntaxVersion']

	def toDictionary(self):
		dictionary = {}
		dictionary['ackResult'] = self.ackResult
		dictionary['ackReason'] = self.ackReason
		dictionary['transferSyntax'] = self.transferSyntax
		dictionary['syntaxVersion'] = self.syntaxVersion
		return dictionary

	def toByteArray(self):
		bytes = bytearray()
		bytes.extend(functions.to16BitLEArray(self.ackResult))
		bytes.extend(functions.to16BitLEArray(self.ackReason))
		bytes.extend(self.transferSyntax.bytes_le)
		bytes.extend(functions.to32BitLEArray(self.syntaxVersion))
		return bytes