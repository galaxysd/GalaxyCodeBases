import binascii
import datetime
import filetimes
import kmsPidGenerator
import struct
import uuid

from structure import Structure

class kmsBase:
	class kmsRequestStruct(Structure):
		commonHdr = ()
		structure = (
			('versionMinor',            '<H'),
			('versionMajor',            '<H'),
			('isClientVm',              '<I'),
			('licenseStatus',           '<I'),
			('graceTime',               '<I'),
			('applicationId',           '16s'),
			('skuId',                   '16s'),
			('kmsCountedId' ,           '16s'),
			('clientMachineId',         '16s'),
			('requiredClientCount',     '<I'),
			('requestTime',             '<Q'),
			('previousClientMachineId', '16s'),
			('machineName',             '128s'),
		)

		def getMachineName(self):
			data = self['machineName'].rsplit('\0\0')[0]
			data += (len(data) & 1 and '\0' or '')
			return data.decode('utf-16le')

		def getLicenseStatus(self):
			return kmsBase.licenseStates[self['licenseStatus']] or "Unknown"

	class kmsResponseStruct(Structure):
		commonHdr = ()
		structure = (
			('versionMinor',         '<H'),
			('versionMajor',         '<H'),
			('epidLen',              '<I=len(kmsEpid)+2'),
			('kmsEpid',              'u'),
			('clientMachineId',      '16s'),
			('responseTime',         '<Q'),
			('currentClientCount',   '<I'),
			('vLActivationInterval', '<I'),
			('vLRenewalInterval',    '<I'),
		)

	appIds = {
		uuid.UUID("55C92734-D682-4D71-983E-D6EC3F16059F") : "Windows",
		uuid.UUID("59A52881-A989-479D-AF46-F275C6370663") : "Office 14 (2010)",
		uuid.UUID("0FF1CE15-A989-479D-AF46-F275C6370663") : "Office 15 (2013)",
	}

	skuIds = {
		uuid.UUID("00091344-1ea4-4f37-b789-01750ba6988c") : "Windows Server 2012 R2 Datacenter",
		uuid.UUID("b3ca044e-a358-4d68-9883-aaa2941aca99") : "Windows Server 2012 R2 Standard",
		uuid.UUID("68531fb9-5511-4989-97be-d11a0f55633f") : "Windows Server 2008 R2 Standard",
		uuid.UUID("81671aaf-79d1-4eb1-b004-8cbbe173afea") : "Windows 8.1 Enterprise",
		uuid.UUID("096ce63d-4fac-48a9-82a9-61ae9e800e5f") : "Windows 8.1 Professional WMC",
		uuid.UUID("c06b6981-d7fd-4a35-b7b4-054742b7af67") : "Windows 8.1 Professional",
		uuid.UUID("fe1c3238-432a-43a1-8e25-97e7d1ef10f3") : "Windows 8.1 Core",
		uuid.UUID("a00018a3-f20f-4632-bf7c-8daa5351c914") : "Windows 8 Professional WMC",
		uuid.UUID("a98bcd6d-5343-4603-8afe-5908e4611112") : "Windows 8 Professional",
		uuid.UUID("ae2ee509-1b34-41c0-acb7-6d4650168915") : "Windows 7 Enterprise",
		uuid.UUID("b92e9980-b9d5-4821-9c94-140f632f6312") : "Windows 7 Professional",
		uuid.UUID("cfd8ff08-c0d7-452b-9f60-ef5c70c32094") : "Windows Vista Enterprise",
		uuid.UUID("4f3d1606-3fea-4c01-be3c-8d671c401e3b") : "Windows Vista Business",
	}

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

	def getUUID(self, data):
		return uuid.UUID(bytes_le=data)
	
	def serverLogic(self, kmsRequest):
		if self.config['debug']:
			print "KMS Request Bytes:", binascii.b2a_hex(str(kmsRequest))
			print "KMS Request:", kmsRequest.dump()

		if self.config['verbose']:
			clientMachineId = self.getUUID(kmsRequest['clientMachineId'])
			applicationId = self.getUUID(kmsRequest['applicationId'])
			skuId = self.getUUID(kmsRequest['skuId'])
			requestDatetime = filetimes.filetime_to_dt(kmsRequest['requestTime'])

			# Try and localize the request time, if pytz is available
			try:
				import timezones
				from pytz import utc
				local_dt = utc.localize(requestDatetime).astimezone(timezones.localtz())
			except ImportError:
				local_dt = requestDatetime

			print "     Machine Name: %s" % kmsRequest.getMachineName()
			print "Client Machine ID: %s" % str(clientMachineId)
			print "   Application ID: %s" % self.appIds.get(applicationId, str(applicationId))
			print "           SKU ID: %s" % self.skuIds.get(skuId, str(skuId))
			print "   Licence Status: %s" % kmsRequest.getLicenseStatus()
			print "     Request Time: %s" % local_dt.strftime('%Y-%m-%d %H:%M:%S %Z (UTC%z)')

		return self.createKmsResponse(kmsRequest)

	def createKmsResponse(self, kmsRequest):
		response = self.kmsResponseStruct()
		response['versionMinor'] = kmsRequest['versionMinor']
		response['versionMajor'] = kmsRequest['versionMajor']

		if not self.config["epid"]:
			response["kmsEpid"] = kmsPidGenerator.epidGenerator(kmsRequest['applicationId'], kmsRequest['versionMajor'], self.config["lcid"]).encode('utf-16le')
		else:
			response["kmsEpid"] = self.config["epid"].encode('utf-16le')
		response['clientMachineId'] = kmsRequest['clientMachineId']
		response['responseTime'] = kmsRequest['requestTime']
		response['currentClientCount'] = self.config["CurrentClientCount"]
		response['vLActivationInterval'] = self.config["VLActivationInterval"]
		response['vLRenewalInterval'] = self.config["VLRenewalInterval"]
		if self.config['verbose']:
			print "      Server ePID: %s" % response["kmsEpid"].decode('utf-16le')
		return response

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
