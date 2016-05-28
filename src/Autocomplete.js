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
