/******/ (function(modules) { // webpackBootstrap
/******/ 	// The module cache
/******/ 	var installedModules = {};

/******/ 	// The require function
/******/ 	function __webpack_require__(moduleId) {

/******/ 		// Check if module is in cache
/******/ 		if(installedModules[moduleId])
/******/ 			return installedModules[moduleId].exports;

/******/ 		// Create a new module (and put it into the cache)
/******/ 		var module = installedModules[moduleId] = {
/******/ 			exports: {},
/******/ 			id: moduleId,
/******/ 			loaded: false
/******/ 		};

/******/ 		// Execute the module function
/******/ 		modules[moduleId].call(module.exports, module, module.exports, __webpack_require__);

/******/ 		// Flag the module as loaded
/******/ 		module.loaded = true;

/******/ 		// Return the exports of the module
/******/ 		return module.exports;
/******/ 	}


/******/ 	// expose the modules object (__webpack_modules__)
/******/ 	__webpack_require__.m = modules;

/******/ 	// expose the module cache
/******/ 	__webpack_require__.c = installedModules;

/******/ 	// __webpack_public_path__
/******/ 	__webpack_require__.p = "";

/******/ 	// Load entry module and return exports
/******/ 	return __webpack_require__(0);
/******/ })
/************************************************************************/
/******/ ([
/* 0 */
/***/ function(module, exports, __webpack_require__) {

	window.TagSelect = __webpack_require__(3);


/***/ },
/* 1 */
/***/ function(module, exports, __webpack_require__) {

	var utils;

	if (true) {
	  utils = {
	    escape: _.escape,
	    invert: _.invert,
	    isEqual: _.isEqual,
	    merge: _.merge,
	    uniq: _.uniq
	  };
	} else {
	  utils = {
	    escape: require('lodash/escape'),
	    escapeRegExp: require('lodash/escapeRegExp'),
	    invert: require('lodash/invert'),
	    isEqual: require('lodash/isEqual'),
	    merge: require('lodash/merge'),
	    uniq: require('lodash/uniq')
	  };
	}

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


/***/ },
/* 2 */
/***/ function(module, exports, __webpack_require__) {

	/**
	 * Autocomplete.js
	 *
	 * @author Graphiq, Inc.
	 */

	var utils = __webpack_require__(1);

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
		this.el.parentNode.insertBefore(this.dd, this.el.nextSibling);

		this.bindEvents();
	};

	/* Operations
	**********************************************/

	Autocomplete.prototype.match = function (input) {
		if (!input) return [];

		var escapedInput = utils.escapeRegExp(input),
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


/***/ },
/* 3 */
/***/ function(module, exports, __webpack_require__) {

	/**
	 * TagSelect.js
	 *
	 * @author Graphiq, Inc.
	 */

	var utils = __webpack_require__(1),
		Autocomplete = __webpack_require__(2);

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


/***/ }
/******/ ]);