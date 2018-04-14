/*
 *****
 ***** ATTENTION: This code has been modified from the original *****
 *****            and adapted for use in PsimagLite
 *****
------------------------- begin original notice -------------------------------

Copyright (C) 2004-2008 René Nyffenegger

This source code is provided 'as-is', without any express or implied
warranty. In no event will the author be held liable for any damages
arising from the use of this software.

Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:

1. The origin of this source code must not be misrepresented; you must not
   claim that you wrote the original source code. If you use this source code
   in a product, an acknowledgment in the product documentation would be
   appreciated but is not required.

2. Altered source versions must be plainly marked as such, and must not be
   misrepresented as being the original source code.

3. This notice may not be removed or altered from any source distribution.

René Nyffenegger rene.nyffenegger@adp-gmbh.ch

------------------------- end original notice -------------------------------

***** Original notice is above
***** ATTENTION: This code has been modified from the original *****
*/

#ifndef PSIBASE64_H
#define PSIBASE64_H

#include "Vector.h"
#include "IoSerializerStub.h"

namespace PsimagLite {

class PsiBase64 {

	static const String base64Chars_;

public:

	class Encode {

	public:

		Encode(const String& str)
		{
			encode_(reinterpret_cast<const unsigned char*>(str.c_str()),str.length());
		}

		Encode(unsigned char const* bytesToEncode, unsigned int inLen)
		{
			encode_(bytesToEncode,inLen);
		}

		const String& operator()() const { return buffer_; }

		void serialize(String label, PsimagLite::IoSerializer& ioSerializer) const
		{
			ioSerializer.write(label, buffer_);
		}

		friend std::ostream& operator<<(std::ostream& os, const Encode& encode)
		{
			os<<"#InputStartsHere\n";
			os<<encode.buffer_;
			os<<"\n#InputEndsHere\n";
			return os;
		}

	private:

		void encode_(unsigned char const* bytesToEncode, unsigned int inLen)
		{
			buffer_ = "";
			int i = 0;
			int j = 0;
			unsigned char charArray3[3];
			unsigned char charArray4[4];

			while (inLen--) {
				charArray3[i++] = *(bytesToEncode++);
				if (i == 3) {
					charArray4[0] = (charArray3[0] & 0xfc) >> 2;
					charArray4[1] = ((charArray3[0] & 0x03)<<4) + ((charArray3[1] & 0xf0)>>4);
					charArray4[2] = ((charArray3[1] & 0x0f)<<2) + ((charArray3[2] & 0xc0)>>6);
					charArray4[3] = charArray3[2] & 0x3f;

					for (i = 0; (i <4) ; i++)
						buffer_ += base64Chars_[charArray4[i]];
					i = 0;
				}
			}

			if (i) {
				for (j = i; j < 3; j++)
					charArray3[j] = '\0';

				charArray4[0] = (charArray3[0] & 0xfc) >> 2;
				charArray4[1] = ((charArray3[0] & 0x03) << 4) + ((charArray3[1] & 0xf0) >> 4);
				charArray4[2] = ((charArray3[1] & 0x0f) << 2) + ((charArray3[2] & 0xc0) >> 6);
				charArray4[3] = charArray3[2] & 0x3f;

				for (j = 0; (j < i + 1); j++)
					buffer_ += base64Chars_[charArray4[j]];

				while ((i++ < 3)) buffer_ += '=';
			}
		}

		String buffer_;
	}; // class Encode

	class Decode {

	public:

		Decode(const String& encodedString)
		{
			buffer_ = "";
			int inLen = encodedString.size();
			int i = 0;
			int j = 0;
			int in_ = 0;
			unsigned char charArray4[4], charArray3[3];

			while (inLen-- && (encodedString[in_] != '=') && isBase64(encodedString[in_])) {
				charArray4[i++] = encodedString[in_]; in_++;
				if (i ==4) {
					for (i = 0; i <4; i++)
						charArray4[i] = base64Chars_.find(charArray4[i]);

					charArray3[0] = (charArray4[0] << 2) + ((charArray4[1] & 0x30) >> 4);
					charArray3[1] = ((charArray4[1] & 0xf) << 4) + ((charArray4[2] & 0x3c) >> 2);
					charArray3[2] = ((charArray4[2] & 0x3) << 6) + charArray4[3];

					for (i = 0; (i < 3); i++)
						buffer_ += charArray3[i];
					i = 0;
				}
			}

			if (i) {
				for (j = i; j <4; j++)
					charArray4[j] = 0;

				for (j = 0; j <4; j++)
					charArray4[j] = base64Chars_.find(charArray4[j]);

				charArray3[0] = (charArray4[0] << 2) + ((charArray4[1] & 0x30) >> 4);
				charArray3[1] = ((charArray4[1] & 0xf) << 4) + ((charArray4[2] & 0x3c) >> 2);
				charArray3[2] = ((charArray4[2] & 0x3) << 6) + charArray4[3];

				for (j = 0; (j < i - 1); j++) buffer_ += charArray3[j];
			}
		}

		const String& operator()() const { return buffer_; }

	private:

		static bool isBase64(unsigned char c) {
			return (isalnum(c) || (c == '+') || (c == '/'));
		}

		String buffer_;
	}; // class Decode
}; // class PsiBase64
}

#endif // PSIBASE64_H

