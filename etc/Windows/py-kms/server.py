import argparse
import binascii
import hashlib
import random
import socket
import SocketServer
import struct
import uuid
import rpcBase, rpcBind, rpcRequest

config = {}

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("ip", nargs="?", action="store", default="0.0.0.0", help="The IP address to listen on. The default is \"0.0.0.0\" (all interfaces).", type=str)
	parser.add_argument("port", nargs="?", action="store", default=1688, help="The network port to listen on. The default is \"1688\".", type=int)
	parser.add_argument("-e", "--epid", dest="epid", default=None, help="Use this flag to manually specify an epid to use. If no epid is specified, a random epid will be generated.", type=str)
	parser.add_argument("-c", "--client-count", dest="CurrentClientCount", default=26, help="Use this flag to specify the current client count. Default is 26. A number >25 is required to enable activation.", type=int)
	parser.add_argument("-a", "--activation-interval", dest="VLActivationInterval", default=120, help="Use this flag to specify the activation interval (in minutes). Default is 120 minutes (2 hours).", type=int)
	parser.add_argument("-r", "--renewal-interval", dest="VLRenewalInterval", default=1440 * 7, help="Use this flag to specify the renewal interval (in minutes). Default is 10080 minutes (7 days).", type=int)
	parser.add_argument("-v", "--verbose", dest="verbose", action="store_const", const=True, default=False, help="Enable this flag to turn on verbose output.")
	parser.add_argument("-d", "--debug", dest="debug", action="store_const", const=True, default=False, help="Enable this flag to turn on debug output. Implies \"-v\".")
	config.update(vars(parser.parse_args()))
	if config['debug']:
		config['verbose'] = True
	server = SocketServer.TCPServer((config['ip'], config['port']), kmsServer)
	server.timeout = 5
	print "TCP server listening at %s on port %d." % (config['ip'],config['port'])
	server.serve_forever()

class kmsServer(SocketServer.BaseRequestHandler):
	def setup(self):
		self.connection = self.request
		print "Connection accepted: %s:%d" % (self.client_address[0],self.client_address[1])

	def handle(self):
		while True:
			# self.request is the TCP socket connected to the client
			try:
				self.data = self.connection.recv(1024)
			except socket.error, e:
				if e[0] == 104:
					break
				else:
					raise
			if self.data == '' or not self.data:
				break
			# self.data = bytearray(self.data.strip())
			# print binascii.b2a_hex(str(self.data))
			localRpcBase = rpcBase.rpcBase(self.data, config)
			packetType = localRpcBase.parseHeader(self.data)['packetType']
			if packetType == localRpcBase.packetType['bindReq']:
				if config['verbose']:
					print "RPC bind request received."
				handler = rpcBind.handler(self.data, config)
			elif packetType == localRpcBase.packetType['request']:
				if config['verbose']:
					print "Received activation request."
				handler = rpcRequest.handler(self.data, config)
			else:
				print "Error: Invalid RPC request type", packetType
				break

			handler.populate()
			res = str(handler.getResponse())
			self.connection.send(res)

			if packetType == localRpcBase.packetType['bindReq']:
				if config['verbose']:
					print "RPC bind acknowledged."
			elif packetType == localRpcBase.packetType['request']:
				if config['verbose']:
					print "Responded to activation request."
				break

	def finish(self):
		self.connection.close()
		print "Connection closed: %s:%d" % (self.client_address[0],self.client_address[1])

if __name__ == "__main__":
	main()
