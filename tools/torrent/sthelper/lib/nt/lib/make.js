var fs             = require('fs');
var path           = require('path');
var OrderedEmitter = require('ordered-emitter');
var b              = require('bncode');
var Buffers        = require('buffers');

var Hasher         = require('./hasher');
var util           = require('./util');
var Torrent        = require('./torrent');
var schema         = require('./schema');
var crcHash = require("crc-hash");
  
// Used to display program info when creating torrents
// and with the cli.
var pkg         = JSON.parse(fs.readFileSync(__dirname + '/../package.json'));
var APP_NAME    = pkg.name;
var APP_VERSION = pkg.version;
var TOLOAD=[],HASH_PERCENT=0,HASH_DONE=0,NOW_HASH='';

/**
 * Makes torrent file object
 *
 * @param {String} announce The announce URL
 * @param {String} dir Directory where the files in the files array are
 * @param {Array.<String>} files Optional
 * @param {Object} options Any options for the torrent go here
 *   {Array.<String>} announceList
 *   {String} comment
 *   {String} name Used in multifile mode. If not present, folder
 *     name will be used.
 *   {Number} pieceLength A power of 2 representing the byte length
 *     of each piece.
 *   {Boolean} private
 *   {Object} moreInfo Typically used to generate a different info hash.
 *   {Number} maxFiles Max files to open at the same time when hashing.
 *   {Number} maxMemory Max amount of memory to allocate when hashing.
 *     Can be bytes or in human readable form such as 300.50MB
 * @param {Function(!Error, Torrent)} callback
 * @return {Hasher}
 */
var make = exports.make = function(announce, dir, files,
                                   options, callback) {
  // First handle all of the optional arguments.
  var typeofFiles = typeof files;
  var typeofOptions = typeof options;

  if (typeofFiles === 'function') {
    callback = files;
    options = {};
    files = ['.'];

  } else if (typeofFiles === 'object' && !Array.isArray(files)) {
    callback = typeofOptions === 'function' ? options : function() {};
    options = files;
    files = ['.'];

  } else {
    if (typeofOptions === 'function') {
      callback = options;
      options = {};

    } else {
      callback = callback || function() {};
      options = options || {};
    }

    if (typeofFiles === 'string') files = [files];
    else if (typeofFiles === 'undefined') files = ['.'];
  }


  var hashOptions = { maxFiles: options.maxFiles, stopOnErr: true };

  // Default piece length is 256kb.
  var pieceLength = options.pieceLength ?
    (1 << options.pieceLength) : 262144;

  // Start hashing files.
  var hasher = new Hasher(dir, files, pieceLength, hashOptions);

  var order  = new OrderedEmitter();
  var buf    = Buffers();
  var info   = {};
  var data, piecesLength, piecesPosition, err;

  hasher.readable = true;

  // Check announce URL.
  if (!util.isURL(announce)) {
    err = new Error('Not a URL: ' + announce);
    process.nextTick(function() {
      hasher.destroy();
      hasher.emit('error', err);
      callback(err);
    });
    return;
  }
  // Check list of files is not empty.
  if (files.length === 0) {
    err = new Error('no files given');
    process.nextTick(function() {
      hasher.destroy();
      hasher.emit('error', err);
      callback(err);
    });
    return;
  }
  // Make metadata object.
  var metadata = {};
  metadata.announce = announce;

  // Check and validate announce list.
  if (options.announceList != null) {
    var msg = schema.announceList(options.announceList);
    if (msg !== null) {
      err = new Error(msg);
      process.nextTick(function() {
        hasher.destroy();
        hasher.emit('error', err);
        callback(err);
      });
      return;
    }
    metadata['announce-list'] = options.announceList;
  }

  // Check comment options.
  if (options.comment != null) metadata.comment = options.comment;

  metadata['created by'] = APP_NAME + ' ' + APP_VERSION;
  metadata['creation date'] = Math.round(Date.now() / 1000);
  metadata.info = info;
	var ofiles=hasher.files;

  // Wait until hasher finishes examining file lengths.
  hasher.on('ready', function() {
  	for(var i=0;i<ofiles.length;i++){
  		TOLOAD.push(ofiles[i])
  	}
  	makecrchash()
    
  });

  // Write last part when hashing is finished.
  hasher.on('end', function() {
  	HASH_DONE++;
	info.pieces = buf.toBuffer();

    // Note that listeners to this `data` event will be called before
    // listeners to the `end` event since its emitted inside this listener.
    hasher.emit('data', data.slice(piecesPosition + piecesLength));
    callback(null, new Torrent(metadata));
    
  });


  // listen for hash events
  hasher.on('hash', function(index, hash) {
    order.emit('hash', { order: index, hash: hash });
  });
  
  hasher.on('progress', function(percent, speed, avg) {
    HASH_PERCENT=percent;
  });

  // Emit hash events to memStream in correct order.
  order.on('hash', function(obj) {
  	
    buf.push(obj.hash);
    hasher.emit('data', obj.hash);
  });


  // Handle stream errors.
  hasher.on('error', function(err) {
    	callback(err);
  });
  make_basicinfo=function (){
  		NOW_HASH='Torrent Pieces'
  		// Multi file mode.
	    if (options.multimode || hasher.files.length > 1) {
	      info.files = ofiles;
	      info.name = options.name || path.basename(dir);

	    // Single file mode.
	    } else {
	      info.length = ofiles[0].length;
	      info.name = path.join.apply(null, ofiles[0].path);
	      info.hash = ofiles[0].hash;
	    }
	    info['piece length'] = pieceLength;
	    info.pieces = null;
	    if (options.private) info.private = 1;

	    // Generate fake pieces to encode them.
	    piecesLength = 20 * hasher.pieces;
	    info.pieces = new Buffer(piecesLength);

	    // Add additional moreInfo to info.
	    for (var moreKey in options.moreInfo) {
	      // Only add moreInfo if it doesn't overwrite info.
	      if (!Object.prototype.hasOwnProperty.call(info, moreKey)) {
	        info[moreKey] = options.moreInfo[moreKey];
	      }
	    }

	    // Bencode data.
	    data = new Buffer(b.encode(metadata), 'binary');
	    var strData = data.toString();
	    var dict = '6:pieces' + piecesLength + ':';
	    piecesPosition = strData.indexOf(dict) + dict.length;
	    piecesPosition = Buffer.byteLength(strData.substr(0, piecesPosition));

	    // Write first part of torrent to stream.
	    process.nextTick(function() {
	      hasher.emit('data', data.slice(0, piecesPosition));

	    });
	    setTimeout(function(){
	    	console.log('making pieces ...')
	    	hasher._start();
	    },1000)
  	}
  	makecrchash=function(){
		var ofile=TOLOAD.shift()
		if(!ofile){
			make_basicinfo()
			return;
		}
		var ofilepath = path.join(dir, path.join.apply(null, ofile.path));
		var rs = fs.createReadStream(ofilepath);
		//console.log('CRC hashing '+ofilepath)
		NOW_HASH=ofilepath
	    var hash = crcHash.createHash("crc32");
	    rs.on('data', hash.update.bind(hash));
	    rs.on('end', function () {
	    	HASH_DONE++;
	        ofile.hash=hash.digest('hex').toUpperCase();
	        //console.log('CRC32 : '+ofile.hash)
	        makecrchash()
	    });
	    rs.on('error', function (e) {
	        ofile.hash='error reading';
	        console.log('error hashing : '+ofilepath)
	        makecrchash()
	    });
	}
	var TIMER=setInterval(function(){
		console.log('Total : '+HASH_DONE+'/'+(ofiles.length+1))
		console.log(NOW_HASH+' : '+HASH_PERCENT+'%');
		
	},1000)
  	return hasher;
};


/**
 * Calls make and creates a write stream with the returned read stream.
 *
 * @param {String} output
 * @param {String} announce
 * @param {String} dir
 * @param {Array.<String>} files
 * @param {Object} options
 * @param {Function(!Error, Torrent)} callback
 * @return {Hasher}
 */
exports.makeWrite = function(output, announce, dir, files, options, callback) {
  var rs = make(announce, dir, files, options, callback);
  rs.pipe(fs.createWriteStream(output));
  return rs;
};


