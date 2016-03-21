import binascii
import datetime
import functions
import kmsPidGenerator
import struct
import uuid

class kmsBase:
	licenseStates = {
		0 : "Unlicensed",
		1 : "Activated",
		2 : "Grace Period",
		3 : "Out-of-Tolerance Grace Period",
		4 : "Non-Genuine Grace Period",
		5 : "Notifications Mode",
		6 : "Extended Grace Period",
	}

	licenseStatesEnum = {
		'unlicensed' : 0,
		'licensed' : 1,
		'oobGrace' : 2,
		'ootGrace' : 3,
		'nonGenuineGrace' : 4,
		'notification' : 5,
		'extendedGrace' : 6
	}

	errorCodes = {
		'SL_E_VL_NOT_WINDOWS_SLP' : 0xC004F035,
		'SL_E_VL_NOT_ENOUGH_COUNT' : 0xC004F038,
		'SL_E_VL_BINDING_SERVICE_NOT_ENABLED' : 0xC004F039,
		'SL_E_VL_INFO_PRODUCT_USER_RIGHT' : 0x4004F040,
		'SL_I_VL_OOB_NO_BINDING_SERVER_REGISTRATION' : 0x4004F041,
		'SL_E_VL_KEY_MANAGEMENT_SERVICE_ID_MISMATCH' : 0xC004F042,
		'SL_E_VL_MACHINE_NOT_BOUND' : 0xC004F056
	}

	unknownBytes = bytearray([ 0x00, 0x00, 0x02, 0x00 ])

	unknownDataSize = 8

	def __init__(self, data, config):
		self.data = data
		self.config = config

	def getConfig(self):
		return self.config

	def getOptions(self):
		return self.config

	def getData(self):
		return self.data

	def getResponse(self):
		return ''

	def getResponsePadding(self, bodyLength):
		if bodyLength % 8 == 0:
			paddingLength = 0
		else:
			paddingLength = 8 - bodyLength % 8
		padding = bytearray(paddingLength)
		return padding

	def serverLogic(self, data):
		kmsRequest = self.parseKmsRequest(data)
		if self.config['verbose']:
			fileTimeInt = struct.unpack("<Q", str(kmsRequest['requestTime']))[0]
			posixTimeInt = (fileTimeInt - 116444736000000000)/10000000
			requestDatetime = datetime.datetime.utcfromtimestamp(posixTimeInt)
			requestDatetimeString = requestDatetime.strftime('%c (UTC)')
			print "Machine Name: %s" % kmsRequest['machineNameString']
			print "CMID: %s" % str(kmsRequest['clientMachineId'])
			print "Application ID: %s" % str(kmsRequest['applicationId'])
			print "SKU ID: %s" % str(kmsRequest['skuId'])
			print "Licence Status: %s" % kmsRequest['licenseStatusString']
			print "Request Time: %s" % requestDatetimeString
		kmsResponse = self.createKmsResponse(kmsRequest)
		epidbuffer = bytearray((kmsResponse['kmsEpid']+'\0').encode('utf-16le'))

		responseArray = bytearray()
		responseArray.extend(functions.to16BitLEArray(kmsRequest['versionMinor']))
		responseArray.extend(functions.to16BitLEArray(kmsRequest['versionMajor']))
		responseArray.extend(functions.to32BitLEArray(len(epidbuffer)))
		responseArray.extend(epidbuffer)
		responseArray.extend(kmsResponse['clientMachineId'].bytes_le)
		responseArray.extend(kmsRequest['requestTime'])
		responseArray.extend(functions.to32BitLEArray(kmsResponse['currentClientCount']))
		responseArray.extend(functions.to32BitLEArray(kmsResponse['vLActivationInterval']))
		responseArray.extend(functions.to32BitLEArray(kmsResponse['vLRenewalInterval']))
		return responseArray

	def createKmsResponse(self, kmsRequest):
		response = {}
		if not self.config["epid"]:
			response["kmsEpid"] = kmsPidGenerator.epidGenerator(kmsRequest['applicationId'], kmsRequest['versionMajor'])
		else:
			response["kmsEpid"] = self.config["epid"]
		response['clientMachineId'] = kmsRequest['clientMachineId']
		response['requestTime'] = kmsRequest['requestTime']
		response['currentClientCount'] = self.config["CurrentClientCount"]
		response['vLActivationInterval'] = self.config["VLActivationInterval"]
		response['vLRenewalInterval'] = self.config["VLRenewalInterval"]
		if self.config['verbose']:
			print "Server ePID: %s" % response["kmsEpid"]
		return response

	def parseKmsRequest(self, data):
		if self.config['debug']:
			print "KMS Request Bytes:", binascii.b2a_hex(str(data))
		kmsRequest = {}
		kmsRequest['versionMinor'] = struct.unpack_from('<H', str(data), 0)[0]
		kmsRequest['versionMajor'] = struct.unpack_from('<H', str(data), 2)[0]
		kmsRequest['isClientVm'] = struct.unpack_from('<I', str(data), 4)[0]
		kmsRequest['licenseStatus'] = struct.unpack_from('<I', str(data), 8)[0]
		kmsRequest['graceTime'] = struct.unpack_from('<I', str(data), 12)[0]
		kmsRequest['applicationId'] = uuid.UUID(bytes_le=str(data[16:16 + 16]))
		kmsRequest['skuId'] = uuid.UUID(bytes_le=str(data[32:32 + 16]))
		kmsRequest['kmsCountedId'] = uuid.UUID(bytes_le=str(data[48:48 + 16]))
		kmsRequest['clientMachineId'] = uuid.UUID(bytes_le=str(data[64:64 + 16]))
		kmsRequest['requiredClientCount'] = struct.unpack_from('<I', str(data), 80)[0]
		kmsRequest['requestTime'] = data[84:84 + 8] # data.ReadUInt64()
		kmsRequest['previousClientMachineId'] = data[92:92 + 16]
		kmsRequest['machineName'] = data[108:108 + 64]

		# translate to human readable
		kmsRequest['machineNameString'] = str(kmsRequest['machineName']).rsplit('\0\0')[0].decode('utf-16le')
		kmsRequest['licenseStatusString'] = self.licenseStates[kmsRequest['licenseStatus']] or "Unknown"

		if self.config['debug']:
			print "KMS Request:", kmsRequest
		return kmsRequest

	def parseVersion(self, data):
		return {
			'versionMajor' : struct.unpack_from('<H', str(data), self.unknownDataSize + 2)[0],
			'versionMinor' : struct.unpack_from('<H', str(data), self.unknownDataSize + 0)[0]
		}

import kmsRequestV4, kmsRequestV5, kmsRequestV6, kmsRequestUnknown

def generateKmsResponseData(data, config):
	localKmsBase = kmsBase(data, config)
	version = localKmsBase.parseVersion(data)['versionMajor']
	currentDate = datetime.datetime.now().ctime()

	if version == 4:
		print "Received V%d request on %s." % (version, currentDate)
		messagehandler = kmsRequestV4.kmsRequestV4(data, config)
		messagehandler.executeRequestLogic()
	elif version == 5:
		print "Received V%d request on %s." % (version, currentDate)
		messagehandler = kmsRequestV5.kmsRequestV5(data, config)
		messagehandler.executeRequestLogic()
	elif version == 6:
		print "Received V%d request on %s." % (version, currentDate)
		messagehandler = kmsRequestV6.kmsRequestV6(data, config)
		messagehandler.executeRequestLogic()
	else:
		print "Unhandled KMS version.", version
		messagehandler = kmsRequestUnknown.kmsRequestUnknown(data, config)
	return messagehandler.getResponse()
