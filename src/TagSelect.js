/**
 * TagSelect.js
 *
 * @author Graphiq, Inc.
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
