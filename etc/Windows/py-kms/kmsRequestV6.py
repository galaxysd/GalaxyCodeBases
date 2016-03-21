import aes
import binascii
import functions
import hashlib
import hmac
import random
import struct
from kmsBase import kmsBase

class kmsRequestV6(kmsBase):
	key = bytearray([ 0xA9, 0x4A, 0x41, 0x95, 0xE2, 0x01, 0x43, 0x2D, 0x9B, 0xCB, 0x46, 0x04, 0x05, 0xD8, 0x4A, 0x21 ])

	def executeRequestLogic(self):
		request = self.parseRequest()
		self.requestData = request

		decrypted = self.decryptRequest(request)
		self.requestData['SaltS'] = decrypted['SaltS'] # Dirty hack part une

		localKmsBase = kmsBase(self.data, self.config)
		responseBuffer = localKmsBase.serverLogic(decrypted['request'])

		encrypted = self.encryptResponse(request, decrypted, responseBuffer)

		data = self.generateResponse(encrypted['response'])

		self.finalResponse = self.generateResponseArray(data)

	def decryptRequest(self, request):
		# SaltS
		SaltS = bytearray(random.getrandbits(8) for i in range(16))

		moo = aes.AESModeOfOperation()
		moo.aes.v6 = True
		decrypted = moo.decrypt(SaltS, 16, moo.modeOfOperation["CBC"], self.key, moo.aes.keySize["SIZE_128"], SaltS)
		#decrypted = aes.strip_PKCS7_padding(decrypted)

		# DSaltS
		DSaltS = functions.bufferToByteArray(decrypted, 0, 16)

		# SaltC
		iv = request['salt']
		SaltC = iv

		encrypted = bytearray()
		encrypted.extend(request['salt'])
		encrypted.extend(request['encryptedRequest'])
		encrypted.extend(request['pad'])

		moo = aes.AESModeOfOperation()
		moo.aes.v6 = True
		decrypted = moo.decrypt(encrypted, 256, moo.modeOfOperation["CBC"], self.key, moo.aes.keySize["SIZE_128"], SaltC)
		decrypted = aes.strip_PKCS7_padding(decrypted)

		# DSaltC
		decryptedSalt = functions.bufferToByteArray(decrypted, 0, 16)
		DSaltC = decryptedSalt

		decryptedRequest = functions.bufferToByteArray(decrypted, 16, len(decrypted) - 16)

		return {
			'salt' : decryptedSalt,
			'request' : decryptedRequest,
			'SaltS' : SaltS,
			'DSaltS' : DSaltS,
			'SaltC' : SaltC,
			'DSaltC' : DSaltC
		}
	
	def encryptResponse(self, request, decrypted, response):
		randomSalt = self.getRandomSalt()
		sha256 = hashlib.sha256()
		sha256.update(str(randomSalt))
		result = sha256.digest()

		randomStuff = bytearray(16)
		for i in range(0,16):
			randomStuff[i] = (decrypted['salt'][i] ^ request['salt'][i] ^ randomSalt[i]) & 0xff

		responsedata = bytearray()
		responsedata.extend(functions.bufferToByteArray(response))
		responsedata.extend(randomStuff)
		responsedata.extend(functions.bufferToByteArray(result))

		# UnknownData
		unknown = bytearray([0x36,0x4F,0x46,0x3A,0x88,0x63,0xD3,0x5F])
		responsedata.extend(unknown)

		# XorSalts
		XorSalts = bytearray(16)
		for i in range (0, 16):
			XorSalts[i] = (decrypted['SaltC'][i] ^ decrypted['DSaltC'][i]) & 0xff
		responsedata.extend(XorSalts)

		# HMacMsg
		HMacMsg = bytearray(230)
		for i in range (0, 16):
			HMacMsg[i] = (decrypted['SaltS'][i] ^ decrypted['DSaltS'][i]) & 0xff
		for i in range (0, 214):
			HMacMsg[i + 16] = responsedata[i]

		# HMacKey
		requestTime = decrypted['request'][-152:-152+8]
		HMacKey = self.getMACKey(requestTime)
		HMac = hmac.new(str(HMacKey), str(HMacMsg), hashlib.sha256)
		digest = HMac.digest()
		responsedata.extend(functions.bufferToByteArray(digest, 16, len(digest) - 16))

		padded = aes.append_PKCS7_padding(responsedata)
		moo = aes.AESModeOfOperation()
		moo.aes.v6 = True
		mode, orig_len, crypted = moo.encrypt(str(padded), moo.modeOfOperation["CBC"], self.key, moo.aes.keySize["SIZE_128"], decrypted['SaltS'])

		encryptedResponse = bytearray(crypted)

		return {
			'response' : encryptedResponse
		}

	def getMACKey(self, timestamp):
		t = struct.unpack("<Q", str(timestamp))[0]
		c1 = 0x00000022816889BD
		c2 = 0x000000208CBAB5ED
		c3 = 0x3156CD5AC628477A

		i1 = (t / c1) & 0xFFFFFFFFFFFFFFFF
		i2 = (i1 * c2) & 0xFFFFFFFFFFFFFFFF
		seed = (i2 + c3) & 0xFFFFFFFFFFFFFFFF

		sha256 = hashlib.sha256()
		sha256.update(struct.pack("<Q", seed))
		digest = sha256.digest()

		return functions.bufferToByteArray(digest, 16, len(digest) - 16)
	
	def parseRequest(self):
		data = self.data
		request = {}
		request['bodyLength1'] = struct.unpack_from('<I', str(data), 0)[0]
		request['bodyLength2'] = struct.unpack_from('<I', str(data), 4)[0]
		request['versionMinor'] = struct.unpack_from('<H', str(data), 8)[0]
		request['versionMajor'] = struct.unpack_from('<H', str(data), 10)[0]
		request['salt'] = functions.bufferToByteArray(data, 12, 16)
		request['encryptedRequest'] = functions.bufferToByteArray(data, 28, len(data) - 8 - 4 - 16 - 4)
		request['pad'] = functions.bufferToByteArray(data, len(data) - 4, 4)
		return request
	
	def getRandomSalt(self):
		return bytearray(random.getrandbits(8) for i in range(16))
	
	def generateResponse(self, encryptedResponse):
		localKmsBase = kmsBase(self.data, self.config)
		bodyLength = 4 + len(self.requestData['SaltS']) + len(encryptedResponse) # Dirty hack part deux
		response = {
			'versionMinor' : self.requestData['versionMinor'],
			'versionMajor' : self.requestData['versionMajor'],
			'bodyLength' : bodyLength,
			'unknown' : localKmsBase.unknownBytes,
			'bodyLength2' : bodyLength,
			'salt' : self.requestData['SaltS'], # Dirty hack part trois
			'encrypted' : encryptedResponse,
			'padding' : localKmsBase.getResponsePadding(bodyLength)
		}

		if self.config['debug']:
			print "KMS V6 Response:", response

		return response
	
	def generateResponseArray(self, data):
		finalResponse = bytearray()
		finalResponse.extend(functions.to32BitLEArray(data['bodyLength']))
		finalResponse.extend(data['unknown'])
		finalResponse.extend(functions.to32BitLEArray(data['bodyLength2']))
		finalResponse.extend(functions.to16BitLEArray(data['versionMinor']))
		finalResponse.extend(functions.to16BitLEArray(data['versionMajor']))
		finalResponse.extend(data['salt'])
		finalResponse.extend(data['encrypted'])
		finalResponse.extend(data['padding'])

		if self.config['debug']:
			print "KMS V6 Response Bytes:", binascii.b2a_hex(str(finalResponse))

		return finalResponse
	
	def getResponse(self):
		return str(self.finalResponse)
