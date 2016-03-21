import binascii
import functions
import kmsBase
import rpcBase
import struct
import uuid

class rpcRequestHandler(rpcBase.rpcBase):
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
		responseArray = bytearray()
		responseArray.append(response['version'])
		responseArray.append(response['versionMinor'])
		responseArray.append(response['packetType'])
		responseArray.append(response['packetFlags'])
		responseArray.extend(functions.to32BitLEArray(response['dataRepresentation']))
		responseArray.extend(functions.to16BitLEArray(response['fragLength']))
		responseArray.extend(functions.to16BitLEArray(response['authLength']))
		responseArray.extend(functions.to32BitLEArray(response['callId']))

		responseArray.extend(functions.to32BitLEArray(response['allocHint']))
		responseArray.append(response['contextId'])
		responseArray.append(response['cancelCount'])
		responseArray.extend(functions.to16BitLEArray(response['opnum']))
		responseArray.extend(response['data'])

		if self.config['debug']:
			print "RPC Message Response Bytes:", binascii.b2a_hex(str(responseArray))

		return responseArray

	def getResponse(self):
		return str(self.responseArray)
