(function e(t,n,r){function s(o,u){if(!n[o]){if(!t[o]){var a=typeof require=="function"&&require;if(!u&&a)return a(o,!0);if(i)return i(o,!0);var f=new Error("Cannot find module '"+o+"'");throw f.code="MODULE_NOT_FOUND",f}var l=n[o]={exports:{}};t[o][0].call(l.exports,function(e){var n=t[o][1][e];return s(n?n:e)},l,l.exports,e,t,n,r)}return n[o].exports}var i=typeof require=="function"&&require;for(var o=0;o<r.length;o++)s(r[o]);return s})({1:[function(require,module,exports){
window.TagSelect = require('./src/TagSelect');

},{"./src/TagSelect":9}],2:[function(require,module,exports){

module.exports = function equal(arr1, arr2) {
  var length = arr1.length
  if (length !== arr2.length) return false
  for (var i = 0; i < length; i++)
    if (arr1[i] !== arr2[i])
      return false
  return true
}

},{}],3:[function(require,module,exports){

'use strict';

module.exports =

    /**
     * Takes in a primitive [keys] and returns a {key:key}
     *  
     * @param array [key1, key2 ... keyn]
     */
    function(arr) {
        return arr.reduce(function (keymirror, key) {
            if(typeof key != 'object' && typeof key != 'function')
                keymirror[key] = key;
            return keymirror;
        }, {});
    };

},{}],4:[function(require,module,exports){
(function (global){
'use strict';

// there's 3 implementations written in increasing order of efficiency

// 1 - no Set type is defined
function uniqNoSet(arr) {
	var ret = [];

	for (var i = 0; i < arr.length; i++) {
		if (ret.indexOf(arr[i]) === -1) {
			ret.push(arr[i]);
		}
	}

	return ret;
}

// 2 - a simple Set type is defined
function uniqSet(arr) {
	var seen = new Set();
	return arr.filter(function (el) {
		if (!seen.has(el)) {
			seen.add(el);
			return true;
		}
	});
}

// 3 - a standard Set type is defined and it has a forEach method
function uniqSetWithForEach(arr) {
	var ret = [];

	(new Set(arr)).forEach(function (el) {
		ret.push(el);
	});

	return ret;
}

// V8 currently has a broken implementation
// https://github.com/joyent/node/issues/8449
function doesForEachActuallyWork() {
	var ret = false;

	(new Set([true])).forEach(function (el) {
		ret = el;
	});

	return ret === true;
}

if ('Set' in global) {
	if (typeof Set.prototype.forEach === 'function' && doesForEachActuallyWork()) {
		module.exports = uniqSetWithForEach;
	} else {
		module.exports = uniqSet;
	}
} else {
	module.exports = uniqNoSet;
}

}).call(this,typeof global !== "undefined" ? global : typeof self !== "undefined" ? self : typeof window !== "undefined" ? window : {})
},{}],5:[function(require,module,exports){
/*!
 * escape-html
 * Copyright(c) 2012-2013 TJ Holowaychuk
 * Copyright(c) 2015 Andreas Lubbe
 * Copyright(c) 2015 Tiancheng "Timothy" Gu
 * MIT Licensed
 */

'use strict';

/**
 * Module variables.
 * @private
 */

var matchHtmlRegExp = /["'&<>]/;

/**
 * Module exports.
 * @public
 */

module.exports = escapeHtml;

/**
 * Escape special characters in the given string of html.
 *
 * @param  {string} string The string to escape for inserting into HTML
 * @return {string}
 * @public
 */

function escapeHtml(string) {
  var str = '' + string;
  var match = matchHtmlRegExp.exec(str);

  if (!match) {
    return str;
  }

  var escape;
  var html = '';
  var index = 0;
  var lastIndex = 0;

  for (index = match.index; index < str.length; index++) {
    switch (str.charCodeAt(index)) {
      case 34: // "
        escape = '&quot;';
        break;
      case 38: // &
        escape = '&amp;';
        break;
      case 39: // '
        escape = '&#39;';
        break;
      case 60: // <
        escape = '&lt;';
        break;
      case 62: // >
        escape = '&gt;';
        break;
      default:
        continue;
    }

    if (lastIndex !== index) {
      html += str.substring(lastIndex, index);
    }

    lastIndex = index + 1;
    html += escape;
  }

  return lastIndex !== index
    ? html + str.substring(lastIndex, index)
    : html;
}

},{}],6:[function(require,module,exports){
'use strict';

var matchOperatorsRe = /[|\\{}()[\]^$+*?.]/g;

module.exports = function (str) {
	if (typeof str !== 'string') {
		throw new TypeError('Expected a string');
	}

	return str.replace(matchOperatorsRe, '\\$&');
};

},{}],7:[function(require,module,exports){
'use strict';
/* eslint-disable no-unused-vars */
var hasOwnProperty = Object.prototype.hasOwnProperty;
var propIsEnumerable = Object.prototype.propertyIsEnumerable;

function toObject(val) {
	if (val === null || val === undefined) {
		throw new TypeError('Object.assign cannot be called with null or undefined');
	}

	return Object(val);
}

function shouldUseNative() {
	try {
		if (!Object.assign) {
			return false;
		}

		// Detect buggy property enumeration order in older V8 versions.

		// https://bugs.chromium.org/p/v8/issues/detail?id=4118
		var test1 = new String('abc');  // eslint-disable-line
		test1[5] = 'de';
		if (Object.getOwnPropertyNames(test1)[0] === '5') {
			return false;
		}

		// https://bugs.chromium.org/p/v8/issues/detail?id=3056
		var test2 = {};
		for (var i = 0; i < 10; i++) {
			test2['_' + String.fromCharCode(i)] = i;
		}
		var order2 = Object.getOwnPropertyNames(test2).map(function (n) {
			return test2[n];
		});
		if (order2.join('') !== '0123456789') {
			return false;
		}

		// https://bugs.chromium.org/p/v8/issues/detail?id=3056
		var test3 = {};
		'abcdefghijklmnopqrst'.split('').forEach(function (letter) {
			test3[letter] = letter;
		});
		if (Object.keys(Object.assign({}, test3)).join('') !==
				'abcdefghijklmnopqrst') {
			return false;
		}

		return true;
	} catch (e) {
		// We don't expect any of the above to throw, but better to be safe.
		return false;
	}
}

module.exports = shouldUseNative() ? Object.assign : function (target, source) {
	var from;
	var to = toObject(target);
	var symbols;

	for (var s = 1; s < arguments.length; s++) {
		from = Object(arguments[s]);

		for (var key in from) {
			if (hasOwnProperty.call(from, key)) {
				to[key] = from[key];
			}
		}

		if (Object.getOwnPropertySymbols) {
			symbols = Object.getOwnPropertySymbols(from);
			for (var i = 0; i < symbols.length; i++) {
				if (propIsEnumerable.call(from, symbols[i])) {
					to[symbols[i]] = from[symbols[i]];
				}
			}
		}
	}

	return to;
};

},{}],8:[function(require,module,exports){
/**
 * Autocomplete.js
 *
 * @author FindTheBest, Inc.
 */

var utils = require('./utils');

/* Constants
**********************************************/

var UP_KEY = 38,
	DOWN_KEY = 40,
	ARROW_KEYS = [UP_KEY, DOWN_KEY],	// Arrow keys
	HIDE_KEYS = [37, 39],				// Keys that should hide suggestions

	TERM_TPL = function (item) {
		var escapedMatch = utils.escape(item.match);
		return (''
		+ '<li class="' + item.ns + '-option' + (item.active ? ' active' : '') + '" '
		+		'data-insert="' + escapedMatch + '">'
		+	escapedMatch
		+ '</li>');
	},

	DEFAULT_SETTINGS = {
		ns: 'tasg',						// CSS/Event namespace
		returnKeys: [13, 9],			// Enter, Tab
		options: [],					// Selectable tokens
		minChars: 0,					// # characters until suggestions appear.
		matchAll: false					// If true, input may also match middle/end of an option
	};

/* Init
**********************************************/

var Autocomplete = function (el, options) {
	this.el = el;
	this.options = utils.merge({}, DEFAULT_SETTINGS, options);
	this.o = this.options;
	this.isVisible = false;
	this.beginOffset = -1;
	this.caretOffset = -1;

	this.dd = document.createElement('ul');
	this.dd.classList.add(this.o.ns + '-dropdown');
	this.el.insertAdjacentElement('afterend', this.dd);

	this.bindEvents();
};

/* Operations
**********************************************/

Autocomplete.prototype.match = function (input) {
	if (!input) return [];

	var escapedInput = utils.escapeStringRegexp(input),
		regStart = new RegExp('^' + escapedInput, 'i'),
		regAll = new RegExp(escapedInput, 'i'),
		matchedStart = [],
		matchedAll = [];

	// Find matches
	for (var word, i = 0; i < this.o.options.length; i++) {
		word = this.o.options[i];
		if (word === input) {
			continue;
		} else if (word.match(regStart)) {
			matchedStart.push(word);
		} else if (this.o.matchAll && word.match(regAll)) {
			matchedAll.push(word);
		}
	}

	return matchedStart.concat(matchedAll);
};

Autocomplete.prototype.select = function (node) {
	if (!this.isVisible) return;
	var inputLength = this.caretOffset - this.beginOffset,
		completion = node.getAttribute('data-insert');

	// Have to delete input and replace it, rather than appending
	//	what's missing. E.g. we might need to replace 'state' with 'State'.
	for (var i = 0; i < inputLength; i++) {
		document.execCommand('delete');
	}
	document.execCommand('insertText', false, completion);
};

Autocomplete.prototype.show = function (matches) {
	matches = matches.map(function (match, i) {
		return { ns: this.o.ns, active: i === 0, match: match };
	}.bind(this));
	this.dd.innerHTML = matches.map(TERM_TPL).join('');
	this.dd.style.left = (this.el.offsetWidth + 1) + 'px';
	this.dd.classList.add('visible');
	this.isVisible = true;
};

Autocomplete.prototype.hide = function () {
	this.dd.classList.remove('visible');
	this.isVisible = false;
};

Autocomplete.prototype.setActive = function (direction) {
	var current = this.dd.querySelector('.active');
	var next = current[direction + 'ElementSibling'];
	if (next) {
		current.classList.remove('active');
		next.classList.add('active');
	}
};

Autocomplete.prototype.setOptions = function (options) {
	this.options = utils.merge(this.options, options);
};

/* Keyboard Events
**********************************************/

Autocomplete.prototype.bindEvents = function () {
	this.el.addEventListener('input', this.onInput.bind(this));
	this.el.addEventListener('keydown', this.onKeydown.bind(this));
	this.el.addEventListener('blur', this.hide.bind(this));
	this.dd.addEventListener('mousedown', function (e) {
		if (e.target.classList.contains(this.o.ns + '-option')) {
			this.select(e.target);
		}
	}.bind(this));
};

Autocomplete.prototype.onInput = function () {
	this.caretOffset = utils.getCaretCharacterOffsetWithin(this.el);
	var text = this.el.textContent;

	// Only make suggestions at the end of a word.
	if (this.caretOffset < text.length && text[this.caretOffset].match(/\w/)) {
		return;
	}

	// Find the beginning of the current word
	this.beginOffset = 0;
	for (var i = this.caretOffset - 1; i >= 0; i--) {
		if (text[i].match(/\W/)) {
			this.beginOffset = i + 1;
			break;
		}
	}

	// Get the partial word and pass it to match()
	var word = text.substr(this.beginOffset, this.caretOffset - this.beginOffset),
		matches = this.match(word);

	if (matches.length) {
		this.show(matches);
	} else {
		this.hide();
	}
};

Autocomplete.prototype.onKeydown = function (e) {
	if (!this.isVisible) {
		return;
	} else if (this.o.returnKeys.indexOf(e.which) >= 0) {
		e.preventDefault();
		e.stopImmediatePropagation();
		this.select(this.dd.querySelector('.active'));
		this.hide();
	} else if (ARROW_KEYS.indexOf(e.which) >= 0) {
		e.preventDefault();
		this.setActive(e.which === UP_KEY ? 'previous' : 'next');
	} else if (HIDE_KEYS.indexOf(e.which) >= 0) {
		this.hide();
	}
};

module.exports = Autocomplete;

},{"./utils":10}],9:[function(require,module,exports){
/**
 * TagSelect.js
 *
 * @author FindTheBest, Inc.
 */

var utils = require('./utils'),
	Autocomplete = require('./Autocomplete');

/****************************************
************* CONSTANTS ****************/

var DEFAULT_ATTRIBUTES = {
		autocorrect: 'off',
		autocapitalize: 'off',
		spellcheck: 'false',
		contenteditable: true
	};

/****************************************
************* TagSelect ****************/

function TagSelect (el, options) {
	el = el.jquery ? el[0] : el;
	options = options || {};
	this.el = el;
	this.autocomplete = null;
	this.selected = options.selected || [];
	this.listeners = { change: [] };

	this.setOptions(options.options);

	el.classList.add('tagselect');
	Object.keys(DEFAULT_ATTRIBUTES).forEach(function (attr) {
		el.setAttribute(attr, DEFAULT_ATTRIBUTES[attr]);
	});

	this.render();
	this.postRender();
}

utils.merge(TagSelect.prototype, {

	render: function () {
		this.shouldValidate = !!Object.keys(this.options).length;
		var html = this.selected
			.map(this.template.bind(this))
			.join('');
		this.el.innerHTML = this.trim(html);
	},

	template: function (field) {
		var escapedField = utils.escape(field),
			invalidClass = this.options[field] || !this.shouldValidate ? '' : ' invalid';
		return (''
		+ '<span data-fd="' + escapedField + '" class="tagselect-fd' + invalidClass + '"> '
		+	escapedField
		+ '</span>'
		+ '<span class="tagselect-fd-separator">, </span>');
	},

	postRender: function () {
		this.autocomplete = new Autocomplete(this.el, {options: Object.keys(this.options)});

		// Must bind autocomplete before this keydown event, since
		// 	the listener will blur the input on [Enter], and
		// 	autocomplete may need to catch that event first.
		this.el.addEventListener('blur', this.update.bind(this));
		this.el.addEventListener('keydown', function (e) {
			// Block mysterious script that keeps blurring the input.
			if (e.which === 91 || e.metaKey) {
				e.stopPropagation();
			} else if (e.which === 13) {
				e.preventDefault();
				this.blur();
			}
		});
	},

	get: function () {
		return this.selected;
	},

	set: function (selected) {
		this.selected = utils.uniq(selected);
		this.render();
	},

	setOptions: function (options) {
		if (Array.isArray(options)) {
			this.options = utils.invert(options);
		} else {
			this.options = options || {};
		}
		if (this.autocomplete) {
			this.autocomplete.setOptions({options: Object.keys(this.options)});
		}
	},

	update: function () {
		var selected = this.el.textContent
				.split(/\W+/)
				.filter(function (s) { return !!s; }),
			changed = !utils.isEqual(selected, this.selected);
		this.set(selected);
		if (changed && this.validate()) {
			this.trigger('change');
		}
	},

	validate: function () {
		if (!Object.keys(this.options).length) return true;

		for (var i = 0, fd; (fd = this.selected[i]); i++) {
			if (!this.options[fd]) {
				return false;
			}
		}
		return true;
	},

	trim: function (text) {
		return text.replace(/\n\s*/gm, '');
	},

	on: function (event, callback) {
		this.listeners[event].push(callback);
	},

	trigger: function (event) {
		this.listeners[event].forEach(function (callback) {
			callback();
		});
	}

});

module.exports = TagSelect;

},{"./Autocomplete":8,"./utils":10}],10:[function(require,module,exports){
var utils;

if ("full" === 'lite') {
  utils = {
    escape: _.escape,
    invert: _.invert,
    isEqual: _.isEqual,
    merge: _.merge,
    uniq: _.uniq
  };
} else {
  utils = {
    escape: require('escape-html'),
    invert: require('array-to-map'),
    isEqual: require('array-equal'),
    merge: require('object-assign'),
    uniq: require('array-uniq')
  };
}

utils.escapeStringRegexp = require('escape-string-regexp');

/**
 * Thanks to Tim Down.
 * http://stackoverflow.com/a/4812022/1314762
 *
 * Linebreaks and some CSS are not handled, so if
 *  that becomes an issue we'll likely need to use
 *  the larger (45kb) rangy-core.js plugin.
 */
utils.getCaretCharacterOffsetWithin = function (element) {
  var caretOffset = 0;
  var doc = element.ownerDocument || element.document;
  var win = doc.defaultView || doc.parentWindow;
  var sel;
  if (typeof win.getSelection !== 'undefined') {
    sel = win.getSelection();
    if (sel.rangeCount > 0) {
      var range = win.getSelection().getRangeAt(0);
      var preCaretRange = range.cloneRange();
      preCaretRange.selectNodeContents(element);
      preCaretRange.setEnd(range.endContainer, range.endOffset);
      caretOffset = preCaretRange.toString().length;
    }
  } else if ( (sel = doc.selection) && sel.type !== 'Control') {
    var textRange = sel.createRange();
    var preCaretTextRange = doc.body.createTextRange();
    preCaretTextRange.moveToElementText(element);
    preCaretTextRange.setEndPoint('EndToEnd', textRange);
    caretOffset = preCaretTextRange.text.length;
  }
  return caretOffset;
};

module.exports = utils;

},{"array-equal":2,"array-to-map":3,"array-uniq":4,"escape-html":5,"escape-string-regexp":6,"object-assign":7}]},{},[1]);
