import binascii
import functions
import rpcBase
import struct
import uuid
from rpcBindRequestCtxItem import *
from rpcBindResponseCtxItem import *

class rpcBindHandler(rpcBase.rpcBase):
	ctxItemsTransferSyntax = []
	ctxItemsTransferSyntax.append(uuid.UUID('8a885d04-1ceb-11c9-9fe8-08002b104860'))
	ctxItemsTransferSyntax.append(uuid.UUID('00000000-0000-0000-0000-000000000000'))
	ctxItemsTransferSyntax.append(uuid.UUID('00000000-0000-0000-0000-000000000000'))

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
		response['fragLength'] = 36 + 3 * 24
		response['authLength'] = request['authLength']
		response['callId'] = request['callId']

		response['maxXmitFrag'] = request['maxXmitFrag']
		response['maxRecvFrag'] = request['maxRecvFrag']
		response['assocGroup'] = 0x1063bf3f

		response['secondaryAddressLength'] = 6
		response['secondaryAddress'] = bytearray(functions.stringPad(self.config['port'], '\0', 6, 'right'))
		response['numberOfResults'] = 3

		response['ctxItems'] = []
		item0 = rpcBindResponseCtxItem({
			'ackResult' : 0,
			'ackReason' : 0,
			'transferSyntax' : self.ctxItemsTransferSyntax[0],
			'syntaxVersion' : 2
		})
		response['ctxItems'].append(item0)
		item1 = rpcBindResponseCtxItem({
			'ackResult' : 2,
			'ackReason' : 2,
			'transferSyntax' : self.ctxItemsTransferSyntax[1],
			'syntaxVersion' : 0
		})
		response['ctxItems'].append(item1)
		item2 = rpcBindResponseCtxItem({
			'ackResult' : 3,
			'ackReason' : 3,
			'transferSyntax' : self.ctxItemsTransferSyntax[2],
			'syntaxVersion' : 0
		})
		response['ctxItems'].append(item2)

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
		responseArray.extend(functions.to32BitLEArray(response['dataRepresentation']))
		responseArray.extend(functions.to16BitLEArray(response['fragLength']))
		responseArray.extend(functions.to16BitLEArray(response['authLength']))
		responseArray.extend(functions.to32BitLEArray(response['callId']))

		responseArray.extend(functions.to16BitLEArray(response['maxXmitFrag']))
		responseArray.extend(functions.to16BitLEArray(response['maxRecvFrag']))
		responseArray.extend(functions.to32BitLEArray(response['assocGroup']))

		responseArray.extend(functions.to16BitLEArray(response['secondaryAddressLength']))
		responseArray.extend(response['secondaryAddress'])
		responseArray.extend(functions.to32BitLEArray(response['numberOfResults']))

		for i in range(0, len(response['ctxItems'])):
			responseArray.extend(response['ctxItems'][i].toByteArray())

		if self.config['debug']:
			print "RPC Bind Response Bytes:", binascii.b2a_hex(str(responseArray))

		return responseArray

	def getResponse(self):
		return str(self.responseArray)
