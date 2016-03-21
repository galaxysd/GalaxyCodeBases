import binascii
import struct
import time
from kmsBase import kmsBase
from structure import Structure

from aes import AES

# v4 AES Key
key = [
	0x05, 0x3D, 0x83, 0x07, 0xF9, 0xE5, 0xF0, 0x88, 0xEB, 0x5E, 0xA6, 0x68, 0x6C, 0xF0, 0x37, 0xC7,
	0x5D, 0x08, 0xBE, 0x94, 0xCD, 0x06, 0x3E, 0x38, 0x41, 0xE2, 0x75, 0x1B, 0x7A, 0xEB, 0xF3, 0x21,
	0xA5, 0xC3, 0x14, 0x49, 0x41, 0x2C, 0xC6, 0x9F, 0xA8, 0x3C, 0xAE, 0xED, 0x8A, 0x51, 0x2B, 0x0B,
	0x8D, 0x7B, 0xF0, 0x74, 0xEE, 0xBC, 0x47, 0x39, 0xD0, 0x2C, 0x85, 0x71, 0x25, 0xBA, 0x21, 0x18,
	0x7B, 0x3A, 0x68, 0xBA, 0x38, 0x58, 0x60, 0x3F, 0xDE, 0xF9, 0x7C, 0xF3, 0x79, 0x74, 0xA6, 0xA0,
	0xC2, 0xB8, 0x94, 0x7A, 0x10, 0x5B, 0xCC, 0x10, 0x4F, 0x94, 0x11, 0x0E, 0xFE, 0xE7, 0xED, 0x08,
	0x8F, 0xA0, 0x71, 0xD1, 0x9B, 0xEC, 0x9D, 0x74, 0x8B, 0x99, 0x56, 0xC9, 0xA3, 0xB4, 0xFD, 0x4B,
	0xCE, 0x40, 0xE2, 0xA8, 0x51, 0x91, 0x41, 0x65, 0x8A, 0xCC, 0xFB, 0xFA, 0xF5, 0x8E, 0x9F, 0x81,
	0x29, 0xF8, 0xDF, 0xA1, 0x7C, 0x98, 0xF5, 0x9B, 0xA6, 0x58, 0xAE, 0x70, 0x2B, 0x31, 0x25, 0x07,
	0x31, 0xC7, 0x00, 0xEB, 0xF4, 0x62, 0x12, 0x1F, 0xFF, 0x0B, 0xFB, 0x43, 0xA5, 0xF3, 0x8D, 0x9E,
	0xC5, 0xDD, 0x79, 0x07, 0x64, 0xA0, 0x7E, 0xEA, 0xEC, 0x25, 0xA6, 0xA6, 0x18, 0x38, 0x8B, 0x71,
	0x4A, 0x68, 0xF1, 0xD1, 0x21, 0x7D, 0x57, 0x3B, 0x45, 0xED, 0x08, 0x9D, 0xA9, 0x4D, 0x8F, 0xD6
]

# Xor Buffer
def xorBuffer(source, offset, destination):
	for i in range(0,16):
		destination[i] ^= source[i + offset]

class kmsRequestV4(kmsBase):
	class RequestV4(Structure):
		commonHdr = ()
		structure = (
			('bodyLength1', '<I'),
			('bodyLength2', '<I'),
			('request',     ':', kmsBase.kmsRequestStruct),
			('hash',        '16s'),
			('padding',     ':'),
		)

	class ResponseV4(Structure):
		commonHdr = ()
		structure = (
			('bodyLength1', '<I=len(response) + len(hash)'),
			('unknown',     '!I=0x00000200'),
			('bodyLength2', '<I=len(response) + len(hash)'),
			('response',    ':', kmsBase.kmsResponseStruct),
			('hash',        '16s'),
			('padding',     ':'),
		)

	def executeRequestLogic(self):
		requestData = self.RequestV4(self.data)

		response = self.serverLogic(requestData['request'])
		hash = self.generateHash(bytearray(str(response)))

		self.responseData = self.generateResponse(response, hash)

		time.sleep(1) # request sent back too quick for Windows 2008 R2, slow it down.

	def generateHash(self, message):
		"""
		The KMS v4 hash is a variant of CMAC-AES-128. There are two key differences:
		* The 'AES' used is modified in particular ways:
		  * The basic algorithm is Rjindael with a conceptual 160bit key and 128bit blocks.
		    This isn't part of the AES standard, but it works the way you'd expect.
		    Accordingly, the algorithm uses 11 rounds and a 192 byte expanded key.
		  * The key doesn't appear to be a normal Rjindael key - which in this case would
		    be a 20 byte key expanded to 192 bytes. It could well just be 192 random
		    bytes.
		  * Odd numbered rounds use a modified operation order where addRoundKey is
		    done second, instead of last.
		* The trailing block is not XORed with a generated subkey, as defined in CMAC.
		  This is probably because the subkey generation algorithm is only defined for
		  situations where block and key size are the same.
		"""
		aes = AES()
		aes.v4 = True

		messageSize = len(message)
		lastBlock = bytearray(16) 
		hashBuffer = bytearray(16)

		# MessageSize / Blocksize
		j = messageSize >> 4

		# Remainding bytes
		k = messageSize & 0xf

		# Hash
		for i in range(0,j):
			xorBuffer(message, i << 4, hashBuffer)
			hashBuffer = aes.encrypt(hashBuffer, key, len(key))

		# Bit Padding
		ii = 0
		for i in range(j << 4, k + (j << 4)):
			lastBlock[ii] = message[i]
			ii += 1
		lastBlock[k] = 0x80

		xorBuffer(lastBlock, 0, hashBuffer)
		hashBuffer = aes.encrypt(hashBuffer, key, len(key))

		return str(bytearray(hashBuffer))

	def generateResponse(self, responseBuffer, hash):
		bodyLength = len(responseBuffer) + len(hash)
		response = self.ResponseV4()
		response['response'] = responseBuffer
		response['hash'] = hash
		response['padding'] = self.getResponsePadding(bodyLength)

		if self.config['debug']:
			print "KMS V4 Response:", response.dump()
			print "KMS V4 Response Bytes:", binascii.b2a_hex(str(response))

		return str(response)

	def getResponse(self):
		return self.responseData

	def generateRequest(self, requestBase):
		hash = str(self.generateHash(bytearray(str(requestBase))))

		bodyLength = len(requestBase) + len(hash)

		request = kmsRequestV4.RequestV4()
		request['bodyLength1'] = bodyLength
		request['bodyLength2'] = bodyLength
		request['request'] = requestBase
		request['hash'] = hash
		request['padding'] = self.getResponsePadding(bodyLength)

		if self.config['debug']:
			print "Request V4 Data:", request.dump()
			print "Request V4:", binascii.b2a_hex(str(request))

		return request
