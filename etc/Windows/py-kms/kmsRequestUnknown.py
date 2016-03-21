import functions
from kmsBase import kmsBase

class kmsRequestUnknown(kmsBase):
	def getResponse(self):
		finalResponse = bytearray();
		finalResponse.extend(functions.to32BitLEArray(0));
		finalResponse.extend(functions.to32BitLEArray(0));
		finalResponse.extend(functions.to32BitLEArray(self.errorCodes['SL_E_VL_KEY_MANAGEMENT_SERVICE_ID_MISMATCH']))
		return str(finalResponse)