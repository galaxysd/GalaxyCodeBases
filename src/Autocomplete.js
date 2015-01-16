/**
 * Autocomplete.js
 *
 * @author FindTheBest, Inc.
 */

/* Constants
**********************************************/

var UP_KEY = 38,
	DOWN_KEY = 40,
	ARROW_KEYS = [UP_KEY, DOWN_KEY],	// Arrow keys
	HIDE_KEYS = [37, 39],				// Keys that should hide suggestions

	TERM_TPL = _.template(''
		+ '<li class="<%= ns %>-option<%= active ? " active" : "" %>" '
		+		'data-insert="<%- match %>">'
		+	'<%- match %>'
		+ '</li>'
	),

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
	this.$el = $(el);
	this.options = _.merge({}, DEFAULT_SETTINGS, options);
	this.o = this.options;
	this.isVisible = false;
	this.beginOffset = -1;
	this.caretOffset = -1;
	this.$dd = $('<ul class="' + this.o.ns + '-dropdown"/>');
	this.$dd.insertAfter(this.$el);

	this.bindEvents();
};

/* Operations
**********************************************/

Autocomplete.prototype.match = function (input) {
	if (!input) return $.Deferred().reject();

	var $matches = $.Deferred(),
		escapedInput = RegExp.quote(input),
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

	// Resolve or reject deferred with matches
	matchedAll = matchedStart.concat(matchedAll);
	if (matchedAll.length) {
		return $matches.resolve(matchedAll);
	}
	return $matches.reject(matchedAll);
};

Autocomplete.prototype.select = function ($node) {
	if (!this.isVisible) return;
	var inputLength = this.caretOffset - this.beginOffset,
		completion = $node.data('insert');

	// Have to delete input and replace it, rather than appending
	//	what's missing. E.g. we might need to replace 'state' with 'State'.
	for (var i = 0; i < inputLength; i++) {
		document.execCommand('delete');
	}
	document.execCommand('insertText', false, completion);
};

Autocomplete.prototype.show = function (matches) {
	matches = _.map(matches, _.bind(function (match, i) {
		return { ns: this.o.ns, active: i === 0, match: match };
	}, this));
	this.$dd
		.html(_.map(matches, TERM_TPL))
		.css({left: this.$el.width() + 1})
		.addClass('visible');
	this.isVisible = true;
};

Autocomplete.prototype.hide = function () {
	this.$dd.removeClass('visible');
	this.isVisible = false;
};

Autocomplete.prototype.setActive = function (direction) {
	var $current = this.$dd.children('.active');
	var $next = $current[direction]();
	if ($next.length) {
		$current.removeClass('active');
		$next.addClass('active');
	}
};

Autocomplete.prototype.setOptions = function (options) {
	this.options = $.extend(this.options, options);
};

/* Keyboard Events
**********************************************/

Autocomplete.prototype.bindEvents = function () {
	this.$el.on({
		input: _.bind(this.onInput, this),
		keydown: _.bind(this.onKeydown, this),
		blur: _.bind(this.hide, this)
	});
	this.$dd.on('mousedown', '.' + this.o.ns + '-option',  _.bind(function (e) {
		this.select($(e.target));
	}, this));
};

Autocomplete.prototype.onInput = function () {
	this.caretOffset = getCaretCharacterOffsetWithin(this.el);
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
	var self = this,
		word = text.substr(this.beginOffset, this.caretOffset - this.beginOffset);
	this.match(word).done(function (matches) {
		self.show(matches);
	}).fail(function () {
		self.hide();
	});
};

Autocomplete.prototype.onKeydown = function (e) {
	if (!this.isVisible) {
		return;
	} else if (_.contains(this.o.returnKeys, e.which)) {
		e.preventDefault();
		e.stopImmediatePropagation();
		this.select(this.$dd.find('.active'));
		this.hide();
	} else if (_.contains(ARROW_KEYS, e.which)) {
		e.preventDefault();
		this.setActive(e.which === UP_KEY ? 'prev' : 'next');
	} else if (_.contains(HIDE_KEYS, e.which)) {
		this.hide();
	}
};

/* Utilities
**********************************************/

/**
 * Thanks to Tim Down.
 * http://stackoverflow.com/a/4812022/1314762
 *
 * Linebreaks and some CSS are not handled, so if
 *	that becomes an issue we'll likely need to use
 *	the larger (45kb) rangy-core.js plugin.
 */
function getCaretCharacterOffsetWithin(element) {
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
}
