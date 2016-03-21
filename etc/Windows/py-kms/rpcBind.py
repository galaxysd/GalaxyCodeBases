import binascii
import rpcBase
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
		bytes.extend(bytearray(struct.pack('<H', self.contextId)))
		bytes.extend(bytearray(struct.pack('<H', self.numTransItems)))
		bytes.extend(self.abstractSyntaxUUID.bytes_le)
		bytes.extend(bytearray(struct.pack('<H', self.abstractSyntaxVersionMajor)))
		bytes.extend(bytearray(struct.pack('<H', self.abstractSyntaxVersionMinor)))
		bytes.extend(self.transferSyntaxUUID.bytes_le)
		bytes.extend(bytearray(struct.pack('<I', self.transferSyntaxVersion)))
		return bytes

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
		bytes.extend(bytearray(struct.pack('<H', self.ackResult)))
		bytes.extend(bytearray(struct.pack('<H', self.ackReason)))
		bytes.extend(self.transferSyntax.bytes_le)
		bytes.extend(bytearray(struct.pack('<I', self.syntaxVersion)))
		return bytes

class handler(rpcBase.rpcBase):
	def parseRequest(self):
		data = self.data
		request = self.parseHeader(data)
		request['maxXmitFrag'] = struct.unpack_from('<H', str(data), 16)[0]
		request['maxRecvFrag'] = struct.unpack_from('<H', str(data), 18)[0]
		request['assocGroup'] = struct.unpack_from('<I', str(data), 20)[0]
		request['numCtxItems'] = struct.unpack_from('<I', str(data), 24)[0]

		request['ctxItems'] = []
		for i in range(0, request['numCtxItems']):
			ctxOffset = (i * rpcBindRequestCtxItem.dataLength) + 28
			localRpcBindRequestCtxItem = rpcBindRequestCtxItem(data[ctxOffset:ctxOffset + rpcBindRequestCtxItem.dataLength])
			request['ctxItems'].append(localRpcBindRequestCtxItem.toDictionary())

		if self.config['debug']:
			print "RPC Bind Request:", request

		return request

	def generateResponse(self):
		response = {}
		request = self.requestData

		response['version'] = request['version']
		response['versionMinor'] = request['versionMinor']
		response['packetType'] = rpcBase.rpcBase.packetType['bindAck']
		response['packetFlags'] = rpcBase.rpcBase.packetFlags['firstFrag'] | rpcBase.rpcBase.packetFlags['lastFrag'] | rpcBase.rpcBase.packetFlags['multiplex']
		response['dataRepresentation'] = request['dataRepresentation']
		response['fragLength'] = 36 + request['numCtxItems'] * 24
		response['authLength'] = request['authLength']
		response['callId'] = request['callId']

		response['maxXmitFrag'] = request['maxXmitFrag']
		response['maxRecvFrag'] = request['maxRecvFrag']
		response['assocGroup'] = 0x1063bf3f

		response['secondaryAddressLength'] = 6
		response['secondaryAddress'] = bytearray(str(self.config['port']).ljust(6, '\0'))
		response['numberOfResults'] = request['numCtxItems']

		preparedResponses = {}
		preparedResponses[self.uuidNDR32] = rpcBindResponseCtxItem({
			'ackResult' : 0,
			'ackReason' : 0,
			'transferSyntax' : self.uuidNDR32,
			'syntaxVersion' : 2
		})
		preparedResponses[self.uuidNDR64] = rpcBindResponseCtxItem({
			'ackResult' : 2,
			'ackReason' : 2,
			'transferSyntax' : self.uuidEmpty,
			'syntaxVersion' : 0
		})
		preparedResponses[self.uuidTime] = rpcBindResponseCtxItem({
			'ackResult' : 3,
			'ackReason' : 3,
			'transferSyntax' : self.uuidEmpty,
			'syntaxVersion' : 0
		})

		response['ctxItems'] = []
		for i in range (0, request['numCtxItems']):
			resp = preparedResponses[request['ctxItems'][i]['transferSyntaxUUID']]
			response['ctxItems'].append(resp)

		if self.config['debug']:
			print "RPC Bind Response:", response

		return response

	def generateResponseArray(self):
		response = self.responseData
		responseArray = bytearray()
		responseArray.append(response['version'])
		responseArray.append(response['versionMinor'])
		responseArray.append(response['packetType'])
		responseArray.append(response['packetFlags'])
		responseArray.extend(bytearray(struct.pack('<I', response['dataRepresentation'])))
		responseArray.extend(bytearray(struct.pack('<H', response['fragLength'])))
		responseArray.extend(bytearray(struct.pack('<H', response['authLength'])))
		responseArray.extend(bytearray(struct.pack('<I', response['callId'])))

		responseArray.extend(bytearray(struct.pack('<H', response['maxXmitFrag'])))
		responseArray.extend(bytearray(struct.pack('<H', response['maxRecvFrag'])))
		responseArray.extend(bytearray(struct.pack('<I', response['assocGroup'])))

		responseArray.extend(bytearray(struct.pack('<H', response['secondaryAddressLength'])))
		responseArray.extend(response['secondaryAddress'])
		responseArray.extend(bytearray(struct.pack('<I', response['numberOfResults'])))

		for i in range(0, len(response['ctxItems'])):
			responseArray.extend(response['ctxItems'][i].toByteArray())

		if self.config['debug']:
			print "RPC Bind Response Bytes:", binascii.b2a_hex(str(responseArray))

		return responseArray

	def getResponse(self):
		return str(self.responseArray)
