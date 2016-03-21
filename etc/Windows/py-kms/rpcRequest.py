import binascii
import kmsBase
import rpcBase
import struct
import uuid

class handler(rpcBase.rpcBase):
	def parseRequest(self):
		data = self.data
		request = self.parseHeader(data)

		request['allocHint'] = struct.unpack_from('<I', str(data), 16)[0]
		request['contextId'] = struct.unpack_from('<H', str(data), 20)[0]
		request['opnum'] = struct.unpack_from('<H', str(data), 22)[0]
		request['data'] = data[24:24 + request['allocHint']]

		return request

	def generateResponse(self):
		response = {}
		request = self.requestData

		response['data'] = kmsBase.generateKmsResponseData(request['data'], self.config)
		envelopeLength = len(response['data'])

		response['version'] = request['version']
		response['versionMinor'] = request['versionMinor']
		response['packetType'] = rpcBase.rpcBase.packetType['response']
		response['packetFlags'] = rpcBase.rpcBase.packetFlags['firstFrag'] | rpcBase.rpcBase.packetFlags['lastFrag']
		response['dataRepresentation'] = request['dataRepresentation']
		response['fragLength'] = envelopeLength + 24
		response['authLength'] = request['authLength']
		response['callId'] = request['callId']

		response['allocHint'] = envelopeLength
		response['contextId'] = request['contextId']
		response['cancelCount'] = 0x00
		response['opnum'] = request['opnum']

		if self.config['debug']:
			print "RPC Message Response:", response

		return response

	def generateResponseArray(self):
		response = self.responseData
		responseArray = str()
		responseArray += struct.pack('B', response['version'])
		responseArray += struct.pack('B', response['versionMinor'])
		responseArray += struct.pack('B', response['packetType'])
		responseArray += struct.pack('B', response['packetFlags'])
		responseArray += struct.pack('<I', response['dataRepresentation'])
		responseArray += struct.pack('<H', response['fragLength'])
		responseArray += struct.pack('<H', response['authLength'])
		responseArray += struct.pack('<I', response['callId'])
		responseArray += struct.pack('<I', response['allocHint'])
		responseArray += struct.pack('<H', response['contextId'])
		responseArray += struct.pack('B', response['cancelCount'])
		responseArray += struct.pack('B', 0) # Reserved
		responseArray += response['data']

		if self.config['debug']:
			print "RPC Message Response Bytes:", binascii.b2a_hex(str(responseArray))

		return responseArray

	def getResponse(self):
		return self.responseArray
