


% syntax: CombineBoxes(combine_type, data_type, data1, data2, data3, data4, options)
%           combine_type      -- 'anelastic', 'pair', 'single'
%           data_type         -- 'displacements', 'phases' or 'test'
%           data1             -- displacement array or phase array if present errors as second column
%           data2             -- time stamp array or amplitude array if present errors as second column
%           data3             -- box X positions 
%           data4             -- box Y positions 
%  options: 'getrid', value   -- remove unwanted row of data if the numner of rows in uneven.
%           'sym_type', value -- 'symmetric' or 'antisymmetric' -- how to combine the rows

test = 0;
%check the combination array is correct.

%test1: single boxes -- output array is 1:5 
test = test + 1;
type = 'single';
fprintf('Test %d: %s boxes\n', test, type)
data1 = [];%1:5;
data2 = [];%1:5;
data3 = [ones(5,1) 3*ones(5,1)];
data4 = ones(5,2);
expected = 1:5;
out = CombineBoxes(type, 'test', data1,data2,data3,data4);
fprintf('run: CombineBoxes(''%s'', ''test'', %s, %s, %s, %s));\n', type, mat2str(data1), mat2str(data2), mat2str(data3), mat2str(data4) )
fprintf('   Expected array\n')
disp(expected)
fprintf('   Returned array\n')
disp(out)
if isequal(out,expected)
    fprintf('Arrays are the same. Test %d passed\n', test)
else
    s = sprintf('Error in script at test %d. the arrays are not the same', test);
    error(s)
end

% test 2. anelastic boxes
test = test + 1;
type = 'anelastic';
data1 = [];%1:5;
data2 = [];%1:5;
data3 = [ones(5,1) 3*ones(5,1)];
data4 = ones(5,2);
expected = [1:4; 2:5];
out = CombineBoxes(type, 'test', data1,data2,data3,data4);
fprintf('\n\nTest %d: %s boxes\n', test, type)
fprintf('run: CombineBoxes(''%s'', ''test'', %s, %s, %s, %s));\n', type, mat2str(data1), mat2str(data2), mat2str(data3), mat2str(data4) )
fprintf('   Expected array\n')
disp(expected)
fprintf('   Returned array\n')
disp(out)
if isequal(out,expected)
    fprintf('Arrays are the same. Test %d passed\n', test)
else
    s = sprintf('Error in script at test %d. the arrays are not the same', test);
    error(s)
end

% test 3. paired boxes
test = test + 1;
fprintf('\n\nTest %d: paired boxes\n', test)
type = 'pair';
data1 = [];%1:5;
data2 = [];%1:5;
data3 = [(1:5)' (2:6)'; (1:5)' (2:6)'];
data4 = [ones(5,2); ones(5,2)*5];
expected = [1:5; 6:10];
fprintf('run: CombineBoxes(''%s'', ''test'', %s, %s, %s, %s));\n', type, mat2str(data1), mat2str(data2), mat2str(data3), mat2str(data4) )
out = CombineBoxes(type, 'test', data1,data2,data3,data4);
fprintf('   Expected array\n')
disp(expected)
fprintf('   Returned array\n')
disp(out)
if isequal(out,expected)
    fprintf('Arrays are the same. Test %d passed\n', test)
else
    s = sprintf('Error in script at test %d. the arrays are not the same', test);
    error(s)
end

% test 4. paired boxes, remove end
test = test + 1;
fprintf('\n\nTest %d: paired boxes remove end\n', test)
type = 'pair';
data1 = [];%1:5;
data2 = [];%1:5;
data3 = [1:5 1:5 1:5; 1:5 1:5 1:5]';
data4 = [1:15; 1:15]';
expected = [1:5; 6:10];
fprintf('run: CombineBoxes(''%s'', ''test'', %s, %s, %s, %s, ''get_rid'', ''end''));\n', type, mat2str(data1), mat2str(data2), mat2str(data3), mat2str(data4) )
out = CombineBoxes(type, 'test', data1,data2,data3,data4, 'get_rid', 'end');
fprintf('   Expected array\n')
disp(expected)
fprintf('   Returned array\n')
disp(out)
if isequal(out,expected)
    fprintf('Arrays are the same. Test %d passed\n', test)
else
    s = sprintf('Error in script at test %d. the arrays are not the same', test);
    error(s)
end

% test 5. paired boxes, remove end
test = test + 1;
fprintf('\n\nTest %d: paired boxes remove middle\n', test)
type = 'pair';
data1 = [];%1:5;
data2 = [];%1:5;
data3 = [1:5 1:5 1:5; 1:5 1:5 1:5]';
data4 = [1:15; 1:15]';
expected = [1:5; 11:15];
fprintf('run: CombineBoxes(''%s'', ''test'', %s, %s, %s, %s, ''get_rid'', ''middle''));\n', type, mat2str(data1), mat2str(data2), mat2str(data3), mat2str(data4) )
out = CombineBoxes(type, 'test', data1,data2,data3,data4, 'get_rid', 'middle');
fprintf('   Expected array\n')
disp(expected)
fprintf('   Returned array\n')
disp(out)
if isequal(out,expected)
    fprintf('Arrays are the same. Test %d passed\n', test)
else
    s = sprintf('Error in script at test %d. the arrays are not the same', test);
    error(s)
end

% test 6. paired boxes, remove 1
test = test + 1;
fprintf('\n\nTest %d: paired boxes remove middle\n', test)
type = 'pair';
data1 = [];%1:5;
data2 = [];%1:5;
data3 = [1:5 1:5 1:5; 1:5 1:5 1:5]';
data4 = [1:15; 1:15]';
expected = [6:10; 11:15];
fprintf('run: CombineBoxes(''%s'', ''test'', %s, %s, %s, %s, ''get_rid'', ''1''));\n', type, mat2str(data1), mat2str(data2), mat2str(data3), mat2str(data4) )
out = CombineBoxes(type, 'test', data1,data2,data3,data4, 'get_rid', '1');
fprintf('   Expected array\n')
disp(expected)
fprintf('   Returned array\n')
disp(out)
if isequal(out,expected)
    fprintf('Arrays are the same. Test %d passed\n', test)
else
    s = sprintf('Error in script at test %d. the arrays are not the same', test);
    error(s)
end


% test 7. paired boxes, remove middle, symmetric
test = test + 1;
fprintf('\n\nTest %d: paired boxes remove middle\n', test)
type = 'pair';
data1 = [];%1:5;
data2 = [];%1:5;
data3 = [1:5 1:5 1:5; 1:5 1:5 1:5]';
data4 = [1:15; 1:15]';
expected = [1:5; 11:15];
fprintf('run: CombineBoxes(''%s'', ''test'', %s, %s, %s, %s, ''get_rid'', ''middle'', ''symm_type'', ''symmetric'');\n', type, mat2str(data1), mat2str(data2), mat2str(data3), mat2str(data4) )
out = CombineBoxes(type, 'test', data1,data2,data3,data4, 'get_rid', 'middle', 'symm_type', 'symmetric');
fprintf('   Expected array\n')
disp(expected)
fprintf('   Returned array\n')
disp(out)
if isequal(out,expected)
    fprintf('Arrays are the same. Test %d passed\n', test)
else
    s = sprintf('Error in script at test %d. the arrays are not the same', test);
    error(s)
end


% test 8. paired boxes, remove middle, asymmetric
test = test + 1;
fprintf('\n\nTest %d: paired boxes remove middle\n', test)
type = 'pair';
data1 = [];%1:5;
data2 = [];%1:5;
data3 = [1:5 1:5 1:5; 1:5 1:5 1:5]';
data4 = [1:15; 1:15]';
expected = [1:5; 15:-1:11];
fprintf('run: CombineBoxes(''%s'', ''test'', %s, %s, %s, %s, ''get_rid'', ''middle'', ''symm_type'', ''asymmetric'');\n', type, mat2str(data1), mat2str(data2), mat2str(data3), mat2str(data4) )
out = CombineBoxes(type, 'test', data1,data2,data3,data4, 'get_rid', 'middle', 'symm_type', 'asymmetric');
fprintf('   Expected array\n')
disp(expected)
fprintf('   Returned array\n')
disp(out)
if isequal(out,expected)
    fprintf('Arrays are the same. Test %d passed\n', test)
else
    s = sprintf('Error in script at test %d. the arrays are not the same', test);
    error(s)
end



%test 9. Test phase combinations.
test = test + 1;
fprintf('\n\nTest %d: paired boxes remove middle, get correct phases\n', test)
type = 'pair';
data1 = [1  1   1 1  1  1 1 1 1 1 1 .5 .5 .5 .5];%1:5;  %phases
data2 = [.5 0 -.1 0 .2 .1 0 0 0 0 .1 .1 .1 .1 .1];%1:5;  %amplitudes
data3 = [1:5 1:5 1:5; 1:5 1:5 1:5]';
data4 = [1:15; 1:15]';
expected1 = [0.4 0 .1 0 .2]; %phases
expected2 = [1 NaN pi+1 NaN 1]; %amplitudes
fprintf('run: CombineBoxes(''%s'', ''phases'', %s, %s, %s, %s, ''get_rid'', ''end'', ''symm_type'', ''symmetric'');\n', type, mat2str(data1), mat2str(data2), mat2str(data3), mat2str(data4) )
[ampout phasout] = CombineBoxes(type, 'phases', data1,data2,data3,data4, 'get_rid', 'end', 'symm_type', 'symmetric');
fprintf('   Expected arrays\n')
fprintf('Phases:   '), disp(expected1)
fprintf('\bAmplitudes'), disp(expected2)
fprintf('   Returned arrays\n')
fprintf('Phases:   '), disp(phasout)
fprintf('\bAmplitudes'), disp(ampout)
if isequaln(phasout,expected1) && isequaln(ampout,expected2)
    fprintf('Arrays are the same. Test %d passed\n', test)
else
    s = sprintf('Error in script at test %d. the arrays are not the same', test);
    error(s)
end


%test 10. Test displacment combinations as well as  reference length and horizontal position arrays.
test = test + 1;
fprintf('\n\nTest %d: paired boxes remove middle, get correct displacements\n', test)
type = 'pair';
data1 = repmat(1:15,150,1); %displacements
data2 = 0:150;  %time stamps
data3 = [1:5 1:5 1:5; 2:6 2:6 2:6]';
data4 = [1*ones(5,2); 5*ones(5,2); 20*ones(5,2)]; data4(:,2) = data4(:,2)+1;
expected = -10*ones(150,5);
expectedlengths = ones(1,5)*19;
expectedhorz = 1.5:1:5.5;
fprintf('run: CombineBoxes(''%s'', ''displacements'', %s, %s, %s, %s, ''get_rid'', ''middle'', ''symm_type'', ''symmetric'');\n', type, mat2str(data1), mat2str(data2), mat2str(data3), mat2str(data4) )
[dispout phasout, lengths, horz_pos] = CombineBoxes(type, 'displacements', data1,data2,data3,data4, 'get_rid', 'middle', 'symm_type', 'symmetric');
fprintf('   Expected arrays: data array size is %d, %d\n', size(expected,1), size(expected,2))
fprintf('    lengths:'), disp(expectedlengths);
fprintf('\b    horizontal position:'), disp(expectedhorz);
% disp(expected)
fprintf('   Returned arrays: data array size is %d, %d\n', size(dispout,1), size(dispout,2))
fprintf('    lengths:'), disp(lengths);
fprintf('\b    horizontal position:'), disp(horz_pos);
% disp(dispout)
if isequaln(dispout,expected) & isequaln(lengths,expectedlengths) & isequaln(horz_pos,expectedhorz)
    fprintf('Arrays are the same. Test %d passed\n', test)
else
    s = sprintf('Error in script at test %d. the arrays are not the same', test);
    error(s)
end


