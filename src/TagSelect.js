(function ( _ ) {

	/****************************************
	************* CONSTANTS ****************/

	var DEFAULT_ATTRIBUTES = {
			autocorrect: 'off',
			autocapitalize: 'off',
			spellcheck: 'false'
		},
		FIELD_TPL = ''
			+	'<% _.forEach(selected, function (fd) { %>'
			+		'<span data-fd="<%- fd %>" class="tagselect-fd<%= options[fd] ? \'\' : \' invalid\' %>">'
			+			'<%- fd %>'
			+		'</span>'
			+		'<span class="tagselect-fd-separator">, </span>'
			+	'<% }) %>';

	/****************************************
	************* TagSelect ****************/

	function TagSelect (options) {
		this.$el = options.$el || options.el;
		this.selected = options.selected || [];
		this.template = _.template(FIELD_TPL);
		this.listeners = { change: [] };

		this.setOptions(options.options);

		this.$el
			.attr(DEFAULT_ATTRIBUTES)
			.prop('contenteditable', true)
			.addClass('tagselect');

		this.render();
		this.postRender();
	}

	$.extend(TagSelect.prototype, {

		render: function () {
			return this.$el.html(
				this.trim(this.template(this))
			);
		},

		postRender: function () {
			this.$el.tagcomplete({options: _.keys(this.options)});

			// Must bind autocomplete before this keydown event, since
			// 	the listener will blur the input on [Enter], and
			// 	autocomplete may need to catch that event first.
			this.$el.on({
				blur: _.bind(this.update, this),
				keydown: function (e) {
					// Block mysterious script that keeps blurring the input.
					if (e.which === 91 || e.metaKey) {
						e.stopPropagation();
					} else if (e.which === 13) {
						e.preventDefault();
						$(this).blur();
					}
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
			if ($.isArray(options)) {
				this.options = _.invert(options);
			} else {
				this.options = options || {};
			}
			if (this.$el.data('tagcomplete')) {
				this.$el.tagcomplete().setOptions({options: _.keys(this.options)});
			}
		},

		update: function () {
			var selected = _.filter(this.$el.text().split(/\W+/)),
				changed = !_.isEqual(selected, this.selected);
			this.set(selected);
			if (changed && this.validate()) {
				this.trigger('change');
			}
		},

		validate: function () {
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
}(window._));