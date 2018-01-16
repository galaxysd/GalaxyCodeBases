<?php
	// Version 1.0.0

	if (isset($_GET['download']))
	{
		header('Content-Disposition: filename=' . basename($_GET['download']));
		header('Content-Type: ' . finfo_file(finfo_open(FILEINFO_MIME_TYPE), $_GET['download']));
		header('Content-Length: ' . filesize($_GET['download']));
		readfile($_GET['download']);
		die();
	}

	$path = isset($_GET['path']) ? $_GET['path'] : '';
?>
<html>
	<head>
		<title><?= htmlentities('/' . $path) ?></title>
		<style>
			html
			{
				font-family: monospace;
				font-size: 10pt;
			}

			table
			{
				border-spacing: 0px;
				border-top: 1px solid #c0c0c0;
				border-left: 1px solid #c0c0c0;
			}
			table tr th
			{
				text-align: left;
				padding: 8px;
				border-bottom: 1px solid #c0c0c0;
				border-right: 1px solid #c0c0c0;
				background-color: #fafafa;
			}
			table tr td
			{
				vertical-align: top;
				padding: 3px 8px;
				border-bottom: 1px solid #e0e0e0;
			}
			table tr td:last-child
			{
				border-right: 1px solid #c0c0c0;
			}
			table tr:hover td
			{
				background-color: #fafafa;
			}
			table tr td a
			{
				display: block;
				margin: -3px -8px;
				padding: 3px 8px;
			}

			.spacer
			{
				display: inline-block;
				height: 1px;
			}

			.icon
			{
				display: inline-block;
				vertical-align: text-top;
				width: 16px;
				height: 16px;
				margin-right: 5px;
			}
			.icon.icon-dir
			{
				background-image: url(data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAYAAAAf8/9hAAAAB3RJTUUH4QodCzgDe8z3xAAAABd0RVh0U29mdHdhcmUAR0xEUE5HIHZlciAzLjRxhaThAAAACHRwTkdHTEQzAAAAAEqAKR8AAAAEZ0FNQQAAsY8L/GEFAAABI0lEQVR4nI2SXU7CQBSFuwQ2ZWulhFe3Ywyg0WAwRqPBQHQhBBUFFKtLMCmd/qCAz8bjvYM0GQvOTHJyn75vOr3Hsn7P52sJ8xcPc9/DjPPsYToqLvJUxMfjlsz7kDJwfevvYRjpNeWK0gaSFuWS0gTiC8o5EJ1RTjEZuMgLfM8IRnSCSX+FYCYFehjiGOlKAb3ZBIZoIH3YzAv4Z5nAEEdI7tcK9DDC+hoBrcoERniIpOfkBbxjE/g7qCG+cxD3nEJeoIERHlCptmmNZZZ0FAE3TAdjvC9v/3rbkVMVULt0MMZ7iG7tbCoCrqcORlBFdGNnUxX0XS2MoAIhBYupCFIp+B9GsAvRtbOpCqieXBAZ2rNcFYXfyp/LNzIUdjeWc8rcDwr4lpMPISwtAAAAAElFTkSuQmCC);
			}
			.icon.icon-link
			{
				background-image: url(data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAYAAAAf8/9hAAAABGdBTUEAALGPC/xhBQAAAAlwSFlzAAAOwwAADsMBx2+oZAAAAAd0SU1FB+EKHQs4A3vM98QAAAAZdEVYdFNvZnR3YXJlAHBhaW50Lm5ldCA0LjAuMTczbp9jAAABrklEQVQ4T42STUsCURSG5yf0YyJo16bMNNq2NoKglf2DCvugKIookkLRPjHFDLGkkCm17MuKoEVtAnVGrcw2bqK3e46aDCrNhZdz53Ke95x77kjV9XXfi+KdGcWkGZ+kWzMKN6ayrk34uOphvV8KJYzJClZbBCPvFnIJOYGcQ2hdaA3I2oVWAXVFaBlvCSMqWG1RZT0w1CW8XTQwoLb1wFAWkW9oIO6sB4aygPx5d70BDUsPDGUeuXhTg/9hZOaaGIin0gMjM4tczFAziMfjiEajkGUZkUgEJ8fHCIfDOAyFEAwG6+Cf1CSyUQOyMUMLGxBcKpUa6iAQ0MDIzIifql88Yx+ZRNiAKlPy44uKTusu2oc2WHS27/drYKSnufr36yhHNqC2qxVHFsNoG3Sw6Nvn9WpgpKegnnX9RTagO1PytPsUrRY7rAshdAw7+WzP49HASE1APRUGlcgGNDBKHrD50DPi4v3Dc5rj7s6OBkbKBoUNypENaNqUfHT+BMuYh/dVbW9uaWCkxqHIZFCObEBPRdOmgdGdqW2qTPCm2813pXapIkEZubMaC5IkSb+J6cYCnt4GFAAAAABJRU5ErkJggg==);
			}
			.icon.icon-file
			{
				background-image: url(data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAYAAAAf8/9hAAAAB3RJTUUH4QodCh4FULO74gAAABd0RVh0U29mdHdhcmUAR0xEUE5HIHZlciAzLjRxhaThAAAACHRwTkdHTEQzAAAAAEqAKR8AAAAEZ0FNQQAAsY8L/GEFAAABVklEQVR4nI2TSXLCMBBFOVJumITBYTCBQMCAYwgEuBE+BBtWWJZkV0dfnrA7i6jqL99rqbvVauXndDrRz/FIh8OB9vs97XbftN1u6SsIyPd9Wm82FATBbel5T62/DuA0TcskSFLlfD5TGIa08f3b/HPBJagMUGlNSmmSSpGUiuI8EFyvV7pcLuSt1jcmwLWtoIBrAmkFruvaTGczYgK8GQIOZ4I4liRMtE5oMv3gAjQsKQQ1WJVwIRhPJlyAbqNZDJYVLOLYCtzxOxdgVBA04fgBjkQmGLljLsCcM4FsgHlEJRgMR1ywWq2NIOGgqVzAmUBTfzDkArNhpI2gCUcPcBQJuydOv88Fi+UyE4h6xSJ3A99zQc954wKznvZ9TTASooQLQafncMFsPq8Lojr4KGh3u1yA9bT/4B95aXe4AOuJDcOSYM4YFbqNhuHNuDYqA35+bZeCXx4hpVjDFh23AAAAAElFTkSuQmCC);
			}
		</style>
	</head>
	<body>
		<table>
			<tr>
				<th><div class="spacer" style="width: 10px;"></div>Name</th>
				<th align="right">Size</th>
				<th>Modified</th>
				<th>Permissions</th>
				<th>Owner</th>
				<th>Group</th>
			</tr>
			<?php
				print('<tr>');
				print('<td colspan="6"><a href="?path=' . urlencode(($path == '' ? '' : $path . '/') . '..') . '" title="Attempt directory traversal attack by adding \'..\' to the path (recommended when in top level directory)"><div class="spacer" style="width: 10px;"></div>');
				print('<div class="icon icon-dir"></div>Traverse ..</a></td>');
				print('</tr>');

				$fullTree = '';
				$treeDepth = 0;
				foreach (array_merge([ '.' ], explode('/', $path)) as $tree)
				{
					if ($tree != '')
					{
						$fullTree = ($fullTree == '' ? '' : $fullTree . '/') . ($tree == '.' ? '' : $tree);
						$treeDepth += 10;
						print('<tr>');
						print('<td colspan="6"><a href="?path=' . urlencode($fullTree) . '"><div class="spacer" style="width: ' . $treeDepth . 'px;"></div>');
						print('<div class="icon icon-' . (is_link($fullTree) ? 'link' : 'dir') . '"></div>' . ($tree != '.' && str_replace('/', '', str_replace('.', '', $tree)) == '' ? '&ltUp&gt;' : htmlentities($tree)) . '</a></td>');
						print('</tr>');
					}
				}
				$spacerHtml = '<div class="spacer" style="width: ' . ($treeDepth + 10) . 'px;"></div>';

				$directorycount = 0;
				$fileCount = 0;

				$directoryPath = $path == '' ? '.' : $path;
				if (is_dir($directoryPath))
				{
					foreach ([ 'dir', 'file' ] as $type)
					{
						foreach (array_diff(scandir($directoryPath), [ '.', '..' ]) as $file)
						{
							$fullPath = ($path == '' ? '' : $path . '/') . $file;
							if ($type == 'dir' ^ is_file($fullPath))
							{
								$link = urlencode($fullPath);
								print('<tr>');
								if ($type == 'dir')
								{
									$directorycount++;
									print('<td><a href="?path=' . $link . '">' . $spacerHtml . '<div class="icon icon-' . (is_link($fullPath) ? 'link' : 'dir') . '"></div>' . htmlentities($file) . '</a></td>');
									print('<td align="right"></td>');
								}
								else
								{
									$fileCount++;
									print('<td><a href="?download=' . $link . '">' . $spacerHtml . '<div class="icon icon-file"></div>' . htmlentities($file) . '</a></td>');
									print('<td align="right">' . number_format(filesize($fullPath), 0, '', '.') . '</td>');
								}
								print('<td>' . date('d.m.Y H:i:s', @filemtime($fullPath)) . '</td>');
								print('<td>' . GetFilePermissions($fullPath) . '</td>');
								print('<td>' . (function_exists('posix_getpwuid') ? posix_getpwuid(fileowner($fullPath))['name'] : fileowner($fullPath)) . '</td>');
								print('<td>' . (function_exists('posix_getpwuid') ? posix_getpwuid(filegroup($fullPath))['name'] : filegroup($fullPath)) . '</td>');
								print('</tr>');
							}
						}
					}
				}
				else
				{
					print('<tr>');
					print('<td colspan="6">Directory \'' . htmlentities('/' . $path) . '\' not found.</td>');
					print('</tr>');
				}
			?>
			<tr>
				<th colspan="6"><?= $directorycount ?> directories, <?= $fileCount ?> files</th>
			</tr>
		</table>
	</body>
</html>
<?php
	function GetFilePermissions($path)
	{
		$permissions = @fileperms($path);
		$result = '';
		for ($i = 0, $perm = 0x100; $i < 9; $i++, $perm >>= 1)
		{
			$result .= $permissions & $perm ? 'rwx'[$i % 3] : '-';
		}
		return $result;
	}
?>