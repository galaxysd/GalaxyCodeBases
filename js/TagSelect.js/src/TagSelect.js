/**
 * TagSelect.js
 *
 * @author FindTheBest, Inc.
 */

/****************************************
************* CONSTANTS ****************/

var DEFAULT_ATTRIBUTES = {
		autocorrect: 'off',
		autocapitalize: 'off',
		spellcheck: 'false',
		contenteditable: true
	},
	FIELD_TPL = ''
		+	'<% _.forEach(selected, function (fd) { %>'
		+		'<span data-fd="<%- fd %>" class="tagselect-fd<%= options[fd] || !shouldValidate ? \'\' : \' invalid\' %>">'
		+			'<%- fd %>'
		+		'</span>'
		+		'<span class="tagselect-fd-separator">, </span>'
		+	'<% }) %>';

/****************************************
************* TagSelect ****************/

function TagSelect (el, options) {
	el = el.jquery ? el[0] : el;
	options = options || {};
	this.el = el;
	this.autocomplete = null;
	this.selected = options.selected || [];
	this.template = _.template(FIELD_TPL);
	this.listeners = { change: [] };

	this.setOptions(options.options);

	el.classList.add('tagselect');
	_.forEach(DEFAULT_ATTRIBUTES, function (value, attr) {
		el.setAttribute(attr, value);
	});

	this.render();
	this.postRender();
}

_.merge(TagSelect.prototype, {

	render: function () {
		this.shouldValidate = !!_.size(this.options);
		this.el.innerHTML = this.trim(this.template(this));
	},

	postRender: function () {
		this.autocomplete = new Autocomplete(this.el, {options: _.keys(this.options)});

		// Must bind autocomplete before this keydown event, since
		// 	the listener will blur the input on [Enter], and
		// 	autocomplete may need to catch that event first.
		this.el.addEventListener('blur', _.bind(this.update, this));
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
		this.selected = _.uniq(selected);
		this.render();
	},

	setOptions: function (options) {
		if (_.isArray(options)) {
			this.options = _.invert(options);
		} else {
			this.options = options || {};
		}
		if (this.autocomplete) {
			this.autocomplete.setOptions({options: _.keys(this.options)});
		}
	},

	update: function () {
		var selected = _.filter(this.el.textContent.split(/\W+/)),
			changed = !_.isEqual(selected, this.selected);
		this.set(selected);
		if (changed && this.validate()) {
			this.trigger('change');
		}
	},

	validate: function () {
		if (!_.size(this.options)) return true;

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
		_.forEach(this.listeners[event], function (callback) {
			callback();
		});
	}

});

// Export TagSelect global
window.TagSelect = TagSelect;
