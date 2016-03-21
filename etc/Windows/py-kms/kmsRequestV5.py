import aes
import binascii
import hashlib
import random
import struct
from kmsBase import kmsBase

class kmsRequestV5(kmsBase):
	key = bytearray([ 0xCD, 0x7E, 0x79, 0x6F, 0x2A, 0xB2, 0x5D, 0xCB, 0x55, 0xFF, 0xC8, 0xEF, 0x83, 0x64, 0xC4, 0x70 ])

	def executeRequestLogic(self):
		request = self.parseRequest()
		self.requestData = request
	
		decrypted = self.decryptRequest(request)
	
		localKmsBase = kmsBase(self.data, self.config)
		responseBuffer = localKmsBase.serverLogic(decrypted['request'])
	
		encrypted = self.encryptResponse(request, decrypted, responseBuffer)
	
		data = self.generateResponse(encrypted['response'])
	
		self.finalResponse = self.generateResponseArray(data)
	
	def decryptRequest(self, request):
		encrypted = bytearray()
		encrypted.extend(request['salt'])
		encrypted.extend(request['encryptedRequest'])
		encrypted.extend(request['pad'])

		iv = request['salt']

		moo = aes.AESModeOfOperation()
		moo.aes.v6 = False
		decrypted = moo.decrypt(encrypted, 256, moo.modeOfOperation["CBC"], self.key, moo.aes.keySize["SIZE_128"], iv)
		decrypted = aes.strip_PKCS7_padding(decrypted)

		decryptedSalt = bytearray(decrypted[:16])

		decryptedRequest = bytearray(decrypted[16:])

		return {
			'salt' : decryptedSalt,
			'request' : decryptedRequest
		}
	
	def encryptResponse(self, request, decrypted, response):
		randomSalt = self.getRandomSalt()
		sha256 = hashlib.sha256()
		sha256.update(str(randomSalt))
		result = sha256.digest()

		randomStuff = bytearray(16)
		for i in range(0,16):
			randomStuff[i] = (decrypted['salt'][i] ^ request['salt'][i] ^ randomSalt[i]) & 0xff

		iv = request['salt']

		responsedata = bytearray()
		responsedata.extend(bytearray(response))
		responsedata.extend(randomStuff)
		responsedata.extend(bytearray(result))

		padded = aes.append_PKCS7_padding(responsedata)
		moo = aes.AESModeOfOperation()
		mode, orig_len, crypted = moo.encrypt(str(padded), moo.modeOfOperation["CBC"], self.key, moo.aes.keySize["SIZE_128"], iv)

		encryptedResponse = bytearray(crypted) #[:192]

		return {
			'response' : encryptedResponse
		}
	
	def parseRequest(self):
		data = self.data
		request = {}
		request['bodyLength1'] = struct.unpack_from('<I', str(data), 0)[0]
		request['bodyLength2'] = struct.unpack_from('<I', str(data), 4)[0]
		request['versionMinor'] = struct.unpack_from('<H', str(data), 8)[0]
		request['versionMajor'] = struct.unpack_from('<H', str(data), 10)[0]
		request['salt'] = bytearray(data[12:28])
		request['encryptedRequest'] = bytearray(data[28:-4])
		request['pad'] = bytearray(data[-4:])
		return request
	
	def getRandomSalt(self):
		return bytearray(random.getrandbits(8) for i in range(16))
	
	def generateResponse(self, encryptedResponse):
		localKmsBase = kmsBase(self.data, self.config)
		bodyLength = 4 + len(self.requestData['salt']) + len(encryptedResponse)
		response = {
			'versionMinor' : self.requestData['versionMinor'],
			'versionMajor' : self.requestData['versionMajor'],
			'bodyLength' : bodyLength,
			'unknown' : localKmsBase.unknownBytes,
			'bodyLength2' : bodyLength,
			'salt' : self.requestData['salt'],
			'encrypted' : encryptedResponse,
			'padding' : localKmsBase.getResponsePadding(bodyLength)
		}

		if self.config['debug']:
			print "KMS V5 Response:", response

		return response
	
	def generateResponseArray(self, data):
		finalResponse = bytearray()
		finalResponse.extend(bytearray(struct.pack('<I', data['bodyLength'])))
		finalResponse.extend(data['unknown'])
		finalResponse.extend(bytearray(struct.pack('<I', data['bodyLength2'])))
		finalResponse.extend(bytearray(struct.pack('<H', data['versionMinor'])))
		finalResponse.extend(bytearray(struct.pack('<H', data['versionMajor'])))
		finalResponse.extend(data['salt'])
		finalResponse.extend(data['encrypted'])
		finalResponse.extend(data['padding'])

		if self.config['debug']:
			print "KMS V5 Response Bytes:", binascii.b2a_hex(str(finalResponse))

		return finalResponse
	
	def getResponse(self):
		return str(self.finalResponse)

