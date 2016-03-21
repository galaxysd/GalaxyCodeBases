import functions
import struct

class rpcBase:
	packetType = {
		'request' : 0,
		'ping' : 1,
		'response' : 2,
		'fault' : 3,
		'working' : 4,
		'nocall' : 5,
		'reject' : 6,
		'ack' : 7,
		'clCancel' : 8,
		'fack' : 9,
		'cancelAck' : 10,
		'bindReq' : 11,
		'bindAck' : 12,
		'bindNak' : 13,
		'alterContext' : 14,
		'alterContextResp' : 15,
		'shutdown' : 17,
		'coCancel' : 18,
		'orphaned' : 19
	}

	packetFlags = {
		'firstFrag' : 1, # 0x01
		'lastFrag' : 2, # 0x02
		'cancelPending' : 4, # 0x04
		'reserved' : 8, # 0x08
		'multiplex' : 16, # 0x10
		'didNotExecute' : 32, # 0x20
		'maybe' : 64, # 0x40
		'objectUuid' : 128 # 0x80
	}

	def __init__(self, data, config):
		self.data = data
		self.config = config

	def populate(self):
		self.requestHeader = self.parseHeader(self.data)
		self.requestData = self.parseRequest()
		self.responseData = self.generateResponse()
		self.responseArray = self.generateResponseArray()
		return self

	def getConfig(self):
		return self.config

	def getOptions(self):
		return self.config

	def getData(self):
		return self.data

	def parseRequest(self):
		return {}

	def getResponse(self):
		return ''

	def parseHeader(self, data):
		header = {}
		header['version'] = struct.unpack_from('B', str(data), 0)[0]
		header['versionMinor'] = struct.unpack_from('B', str(data), 1)[0]
		header['packetType'] = struct.unpack_from('B', str(data), 2)[0]
		header['packetFlags'] = struct.unpack_from('B', str(data), 3)[0]
		header['dataRepresentation'] = struct.unpack_from('<I', str(data), 4)[0]
		header['fragLength'] = struct.unpack_from('<H', str(data), 8)[0]
		header['authLength'] = struct.unpack_from('<H', str(data), 10)[0]
		header['callId'] = struct.unpack_from('<I', str(data), 12)[0]
		return header
