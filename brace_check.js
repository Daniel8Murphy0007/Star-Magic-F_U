const fs = require('fs');
const content = fs.readFileSync('index.js', 'utf8');
let braceCount = 0;
let inString = false;
let stringChar = '';
let inComment = false;
let inMultiComment = false;

for (let i = 0; i < content.length; i++) {
    const char = content[i];
    const nextChar = content[i + 1] || '';
    
    if (!inString && !inComment && !inMultiComment) {
        if (char === '/' && nextChar === '/') {
            inComment = true;
            continue;
        }
        if (char === '/' && nextChar === '*') {
            inMultiComment = true;
            i++;
            continue;
        }
    }
    
    if (inComment && char === '\n') {
        inComment = false;
        continue;
    }
    
    if (inMultiComment && char === '*' && nextChar === '/') {
        inMultiComment = false;
        i++;
        continue;
    }
    
    if (inComment || inMultiComment) continue;
    
    if (!inString && (char === '\"' || char === '\'')) {
        inString = true;
        stringChar = char;
        continue;
    }
    
    if (inString && char === stringChar && content[i - 1] !== '\\\\') {
        inString = false;
        stringChar = '';
        continue;
    }
    
    if (inString) continue;
    
    if (char === '{') {
        braceCount++;
    } else if (char === '}') {
        braceCount--;
        if (braceCount < 0) {
            console.log('Negative at position', i);
            break;
        }
    }
}

console.log('Final brace count:', braceCount);
