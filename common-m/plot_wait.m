function cmd = plot_wait(name)
% Present replotting output options and await command
%
% Note that this function works in Octave version 2.0.16, but not in
% Octave version 2.1.57.
% Partial work around: After using this function, do "gset term fig"
% (fig big may not work at all in the later Octave version).
	rotx = 30;
	rotz = 60;
	cmd = 'dummy';
	while (cmd!='')
		cmd = input(['Press a key to continue [f/e/p=save ',name,'.fig/eps/ps, P=print, n=name, x=re-X, y=re-Y, r=rotate-Z, R=rotate-X]...'],'s');
		if cmd=='n'
			name = input('Change file name to: ','s');
		end
		if cmd=='f'
			subcmd = input('Figure plot [n]ormal, [b]ig or [p]ortrait? ','s');
			if (subcmd=='p')
			  __gnuplot_set__ term fig big portrait;
			else
			  if (subcmd=='b')
			    __gnuplot_set__ term fig big;
			  else
			    __gnuplot_set__ term fig;
			  endif
			endif
			eval(['__gnuplot_set__ output "',name,'.fig"']);
			replot;
			__gnuplot_set__ term x11;
			printf('%s %s.%s\n','Saved',name,'fig');
		end
		if cmd=='e'
			__gnuplot_set__ term post eps;
			eval(['__gnuplot_set__ output "',name,'.eps"']);
			replot;
			__gnuplot_set__ term x11;
			printf('%s %s.%s\n','Saved',name,'eps');
		end
		if cmd=='p'
			__gnuplot_set__ size 0.721,0.721;
			__gnuplot_set__ term post portrait;
			eval(['__gnuplot_set__ output "',name,'.ps"']);
			replot;
			__gnuplot_set__ size 1.0,1.0;
			__gnuplot_set__ term x11;
			printf('%s %s.%s\n','Saved',name,'ps');
		end
		if cmd=='P'
			[u_status, u_msg] = unlink(['print-tmp-',name,'.ps']);
			%__gnuplot_set__ size 0.721,0.721;
			__gnuplot_set__ term post landscape;
			eval(['__gnuplot_set__ output "print-tmp-',name,'.ps"']);
			replot;
			__gnuplot_set__ size 1.0,1.0;
			__gnuplot_set__ term x11;
			system(['lp -dhp "print-tmp-',name,'.ps"']);
			printf('Sending print-tmp-%s.ps to spool: ',name);
			for spooldelay=0:5,
				printf('+'); fflush(1);
				sleep(1);
			endfor
			printf('\nPrinting print-tmp-%s.ps\n',name);
			[u_status, u_msg] = unlink(['print-tmp-',name,'.ps']);
		end
		if cmd=='x'
			xrg = input('Set X-range to [default=[-1.0 1.0]]: ');
			if isempty(xrg)
				xrg = [-1.0 1.0];
			endif
			if (length(xrg)!=2)
				xrg = [-1.0 xrg(1)];
			endif
			eval(['__gnuplot_set__ xrange [',num2str(xrg(1)),':',num2str(xrg(2)),'];']);
			replot;
		endif
		if cmd=='y'
			yrg = input('Set Y-range to [default=[-1.0 1.0]]: ');
			if isempty(yrg)
				yrg = [-1.0 1.0];
			endif
			if (length(yrg)!=2)
				yrg = [-1.0 yrg(1)];
			endif
			eval(['__gnuplot_set__ yrange [',num2str(yrg(1)),':',num2str(yrg(2)),'];']);
			replot;
		endif
		if cmd=='r'
			rotz = rotz + 30;
			if (rotz>=360)
				rotz = 0;
			endif
			eval(['__gnuplot_set__ view ',num2str(rotx),', ',num2str(rotz),', 1, 1']);
			replot
		endif
		if cmd=='R'
			rotx = rotx + 30;
			if (rotx>=180)
				rotx = 0;
			endif
			eval(['__gnuplot_set__ view ',num2str(rotx),', ',num2str(rotz),', 1, 1']);
			replot
		endif
	endwhile
	__gnuplot_set__ term x11;
endfunction
