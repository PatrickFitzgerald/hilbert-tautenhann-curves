function packageFolder = getPackagePath()
	thisFilePath = mfilename('fullpath');
	[packageFolder,~,~] = fileparts(thisFilePath);
end